import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.shared.utils import copy_dict
from colabdesign.shared.model import soft_seq
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.af.alphafold.model import model, config

############################################################################
# AF_INPUTS - functions for modifying inputs before passing to alphafold
############################################################################
class _af_inputs:

  def set_seq(self, seq=None, mode=None, **kwargs):
    if self._args["optimize_seq"]:
      self._set_seq(seq=seq, mode=mode, **kwargs)
    
    else:
      seq,_ = self._set_seq(seq=seq, mode=mode, return_values=True, **kwargs)
      if self.protocol == "binder":
        # concatenate target and binder sequence
        seq_target = self._inputs["batch"]["aatype"][:self._target_len]
        seq_target_oh = np_one_hot(seq_target, self._args["alphabet_size"])
        seq_target_oh = np.broadcast_to(seq_target_oh,(self._num, *seq_target_oh.shape))
        seq = np.concatenate([seq_target_oh,seq],1)
      
      if self.protocol in ["fixbb","hallucination","partial"] and self._args["copies"] > 1:
        seq = expand_copies(seq, self._args["copies"], self._args["block_diag"], use_jax=False)

      update_seq(seq, self._inputs, use_jax=False)
      update_aatype(seq, self._inputs, use_jax=False)

  def _update_seq(self, params, inputs, aux, key):
    '''get sequence features'''
    opt = inputs["opt"]
    k1, k2 = jax.random.split(key)
    seq = soft_seq(params["seq"], inputs["bias"], opt, k1, num_seq=self._num,
                   shuffle_first=self._args["shuffle_first"])
    seq = self._fix_pos(seq)
    aux.update({"seq":seq, "seq_pseudo":seq["pseudo"]})
    
    # protocol specific modifications to seq features
    if self.protocol == "binder":
      # concatenate target and binder sequence
      seq_target = jax.nn.one_hot(inputs["batch"]["aatype"][:self._target_len],self._args["alphabet_size"])
      seq_target = jnp.broadcast_to(seq_target,(self._num, *seq_target.shape))
      seq = jax.tree_map(lambda x:jnp.concatenate([seq_target,x],1), seq)
      
    if self.protocol in ["fixbb","hallucination","partial"] and self._args["copies"] > 1:
      seq = jax.tree_map(lambda x:expand_copies(x, self._args["copies"], self._args["block_diag"]), seq)
              
    # update sequence features      
    update_seq(seq, inputs)
    
    # update amino acid sidechain identity
    update_aatype(seq, inputs) 
    inputs["seq"] = aux["seq"]
    
    return seq

  def _fix_pos(self, seq, return_p=False):
    if "fix_pos" in self.opt:
      if "pos" in self.opt:
        seq_ref = jax.nn.one_hot(self._wt_aatype_sub,self._args["alphabet_size"])
        p = self.opt["pos"][self.opt["fix_pos"]]
        fix_seq = lambda x: x.at[...,p,:].set(seq_ref)
      else:
        seq_ref = jax.nn.one_hot(self._wt_aatype,self._args["alphabet_size"])
        p = self.opt["fix_pos"]
        fix_seq = lambda x: x.at[...,p,:].set(seq_ref[...,p,:])
      seq = jax.tree_map(fix_seq, seq)
      if return_p: return seq, p
    return seq

  def _update_template(self, inputs):
    ''''dynamically update template features''' 

    # gather features
    opt, batch = inputs["opt"], inputs["batch"]
    if batch["aatype"].ndim == 1: batch = jax.tree_map(lambda x:x[None], batch)
    T,L = batch["aatype"].shape

    # decide which position to remove sequence and/or sidechains
    rm     = jnp.broadcast_to(inputs.get("rm_template",False),L)
    rm_seq = jnp.where(rm,True,jnp.broadcast_to(inputs.get("rm_template_seq",True),L))
    rm_sc  = jnp.where(rm_seq,True,jnp.broadcast_to(inputs.get("rm_template_sc",True),L))

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
    
    if self.protocol == "partial":
      # figure out what positions to use
      pos = opt["pos"]
      if self._args["repeat"] or self._args["homooligomer"]:
        # if repeat/homooligomer propagate the positions across copies 
        C = self._args["copies"]
        pos = (jnp.repeat(pos,C).reshape(-1,C) + jnp.arange(C) * self._len).T.flatten()
      
      # inject template features
      for k,v in template_feats.items():
        if k in ["template_dgram"]:
          inputs[k] = inputs[k].at[:T,pos[:,None],pos[None,:]].set(v)
        else:
          inputs[k] = inputs[k].at[:T,pos].set(v)

      # remove sidechains (mask anything beyond CB)
      k = "template_all_atom_mask"
      inputs[k] = inputs[k].at[:T,pos,5:].set(jnp.where(rm_sc[:,None],0,inputs[k][:T,pos,5:]))
      inputs[k] = inputs[k].at[:T,pos].set(jnp.where(rm[:,None],0,inputs[k][:T,pos]))

    else:
      # inject template features
      for k,v in template_feats.items():
        inputs[k] = inputs[k].at[:T].set(v)      
      # remove sidechains (mask anything beyond CB)
      k = "template_all_atom_mask"
      inputs[k] = inputs[k].at[:T,:,5:].set(jnp.where(rm_sc[:,None],0,inputs[k][:T,:,5:]))
      inputs[k] = jnp.where(rm[:,None],0,inputs[k])

def update_seq(seq, inputs, use_jax=True):
  '''update the sequence features'''  

  _np = jnp if use_jax else np
  msa = seq["pseudo"] if isinstance(seq, dict) else seq
    
  target_feat = msa[0,:,:20]
  msa = _np.pad(msa,[[0,0],[0,0],[0,22-msa.shape[-1]]])
  inputs.update({"msa":msa, "target_feat":target_feat})

def update_aatype(seq, inputs, use_jax=True):
  _np = jnp if use_jax else np
  if isinstance(seq,dict):
    aatype = seq["pseudo"][0].argmax(-1)
  else:
    aatype = seq[0].argmax(-1)
  inputs["aatype"] = aatype

def np_one_hot(x, alphabet):
  return np.pad(np.eye(alphabet),[[0,1],[0,0]])[x]

def expand_copies(x, copies, block_diag=True, use_jax=True):
  '''
  given msa (N,L,20) expand to (1+N*copies,L*copies,22) if block_diag else (N,L*copies,22)
  '''
  _np = jnp if use_jax else np
  one_hot = jax.nn.one_hot if use_jax else np_one_hot
  if x.shape[-1] < 22:
    x = _np.pad(x,[[0,0],[0,0],[0,22-x.shape[-1]]])
  x = _np.tile(x,[1,copies,1])
  if copies > 1 and block_diag:
    L = x.shape[1]
    sub_L = L // copies
    y = x.reshape((-1,1,copies,sub_L,22))
    block_diag_mask = _np.expand_dims(_np.eye(copies),(0,3,4))
    seq = block_diag_mask * y
    gap_seq = (1-block_diag_mask) * one_hot(_np.repeat(21,sub_L),22)
    y = (seq + gap_seq).swapaxes(0,1).reshape(-1,L,22)
    return _np.concatenate([x[:1],y],0)
  else:
    return x