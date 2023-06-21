import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.shared.utils import copy_dict
from colabdesign.shared.model import soft_seq
from colabdesign.shared.parsers import parse_a3m
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.af.alphafold.model import model, config

############################################################################
# AF_INPUTS - functions for modifying inputs before passing to alphafold
############################################################################
class _af_inputs:

  def set_seq(self, seq=None, mode=None, **kwargs):
    assert self._args["optimize_seq"] == True
    self._set_seq(seq=seq, mode=mode, **kwargs)

  def set_msa(self, msa=None, deletion_matrix=None, a3m_filename=None):
    ''' set msa '''
    assert self._args["optimize_seq"] == False

    if a3m_filename is not None:
      msa, deletion_matrix = parse_a3m(a3m_filename=a3m_filename)

    if msa is None:
      msa = np.zeros((self._num, self._len),int)

    if msa.ndim == 1:
      msa = msa[None]

    if deletion_matrix is None:
      deletion_matrix = np.zeros(msa.shape)

    if self.protocol == "binder":
      # add target sequence
      msa = np.pad(msa,[[0,0],[self._target_len,0]])
      deletion_matrix = np.pad(deletion_matrix,[[0,0],[self._target_len,0]])
      msa[0,:self._target_len] = self._inputs["batch"]["aatype"][:self._target_len]
      msa[1:,:self._target_len,-1] = 1
    
    elif self._args["copies"] > 1:
      msa, deletion_matrix = expand_copies(msa, copies=self._args["copies"], 
        block_diag=self._args["block_diag"], d=deletion_matrix, use_jax=False)

    self._inputs.update({
      "msa":msa,
      "target_feat":np_one_hot(msa[0],20),
      "aatype":msa[0],
      "deletion_matrix":deletion_matrix,
      "msa_mask":np.ones_like(deletion_matrix)
    })
  
  def _update_seq(self, params, inputs, aux, key):
    '''get sequence features'''

    def _fix_pos(seq):
      if "fix_pos" in self.opt:
        seq_ref = jax.nn.one_hot(self._wt_aatype,self._args["alphabet_size"])
        p = self.opt["fix_pos"]
        fix_seq = lambda x: x.at[...,p,:].set(seq_ref[...,p,:])
        seq = jax.tree_map(fix_seq, seq)
      return seq

    opt = inputs["opt"]
    k1, k2 = jax.random.split(key)
    seq = soft_seq(params["seq"], inputs["bias"], opt, k1, num_seq=self._num,
                   shuffle_first=self._args["shuffle_first"])
    seq = _fix_pos(seq)
    aux.update({"seq":seq, "seq_pseudo":seq["pseudo"]})
    
    # protocol specific modifications to seq features
    if self.protocol == "binder":
      # concatenate target and binder sequence
      seq_target = jax.nn.one_hot(inputs["batch"]["aatype"][:self._target_len],self._args["alphabet_size"])
      seq_target = jnp.broadcast_to(seq_target,(self._num, *seq_target.shape))
      seq = jax.tree_map(lambda x:jnp.concatenate([seq_target,x],1), seq)
      
    elif self._args["copies"] > 1:
      seq = jax.tree_map(lambda x:expand_copies(x, self._args["copies"], self._args["block_diag"]), seq)
              
    # update inputs
    msa = seq["pseudo"]
    msa = jnp.pad(msa,[[0,0],[0,0],[0,22-msa.shape[-1]]])      
    inputs.update({
      "msa":msa,
      "target_feat":msa[0,:,:20],
      "aatype":msa[0].argmax(-1),
      "deletion_matrix":jnp.zeros(msa.shape[:2]),
      "msa_mask":jnp.ones(msa.shape[:2])
    })
    
    return seq

  def _update_template(self, inputs):
    ''''dynamically update template features''' 

    # gather features
    opt, batch = inputs["opt"], inputs["batch"]
    if batch["aatype"].ndim == 1: batch = jax.tree_map(lambda x:x[None], batch)
    T,L = batch["aatype"].shape

    # decide which position to remove sequence and/or sidechains
    opt_T = opt["template"]
    rm     = jnp.broadcast_to(inputs.get("rm_template",opt_T["rm"]),L)
    rm_seq = jnp.where(rm,True,jnp.broadcast_to(inputs.get("rm_template_seq",opt_T["rm_seq"]),L))
    rm_sc  = jnp.where(rm_seq,True,jnp.broadcast_to(inputs.get("rm_template_sc",opt_T["rm_sc"]),L))

    # define template features
    template_feats = {"template_aatype":jnp.where(rm_seq,21,batch["aatype"])}

    if "dgram" in batch:
      # use dgram from batch if provided
      template_feats.update({"template_dgram":batch["dgram"]})
      nT,nL = inputs["template_aatype"].shape
      inputs["template_dgram"] = jnp.zeros((nT,nL,nL,39))
      
    if "all_atom_positions" in batch:
      # use coordinates from batch if provided
      template_feats.update({"template_all_atom_positions": batch["all_atom_positions"],
                             "template_all_atom_mask":      batch["all_atom_mask"]})

    # enable templates
    inputs["template_mask"] = inputs["template_mask"].at[:T].set(1)    
    
    # inject template features
    for k,v in template_feats.items():
      inputs[k] = inputs[k].at[:T].set(v)      
    # remove sidechains (mask anything beyond CB)
    k = "template_all_atom_mask"
    inputs[k] = inputs[k].at[:T,:,5:].set(jnp.where(rm_sc[:,None],0,inputs[k][:T,:,5:]))
    inputs[k] = jnp.where(rm[:,None],0,inputs[k])

def np_one_hot(x, alphabet):
  return np.pad(np.eye(alphabet),[[0,1],[0,0]])[x]

def expand_copies(x, copies=1, block_diag=True, d=None, use_jax=True):
  '''
  given msa (N,L,20) expand to (1+N*copies,L*copies,22) if block_diag else (N,L*copies,22)
  '''
  _np = jnp if use_jax else np

  N,L = x.shape[:2]
  if x.ndim == 3:
    x = _np.pad(x,[[0,0],[0,0],[0,22-x.shape[-1]]])
    x = _np.tile(x,[1,copies,1])
  else:
    x = _np.tile(x,[1,copies])
  
  if d is not None:
    d = _np.tile(d,[1,copies])

  if copies > 1 and block_diag:
    i = _np.arange(copies)    
    if x.ndim == 3:
      x_ = x.reshape((N,copies,L,22))
      x_diag = _np.zeros((N,copies,copies,L,22))
      if use_jax:
        x_diag = x_diag.at[...,-1].set(1).at[:,i,i].set(x_)
      else:
        x_diag[...,-1] = 1
        x_diag[:,i,i] = x_
      x_diag = x_diag.swapaxes(0,1).reshape((N*copies,copies*L,22))
    
    else:
      x_ = x.reshape((N,copies,L))
      x_diag = _np.full((N,copies,copies,L),21)
      if use_jax:
        x_diag = x_diag.at[:,i,i].set(x_)
      else:
        x_diag[:,i,i] = x_
      x_diag = x_diag.swapaxes(0,1).reshape((N*copies,copies*L))    
    
    if d is not None:
      d_ = d.reshape((N,copies,L))
      d_diag = _np.zeros((N,copies,copies,L))
      if use_jax:
        d_diag = d_diag.at[:,i,i].set(d_)
      else:
        d_diag[:,i,i] = d_
      d_diag = d_diag.swapaxes(0,1).reshape((N*copies,copies*L))
      d = _np.concatenate([d[:1],d_diag],0)
    
    x = _np.concatenate([x[:1],x_diag],0)
  
  return x if d is None else (x,d)