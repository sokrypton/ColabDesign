import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.af.alphafold.common import residue_constants
from colabdesign.af.alphafold.model import model

############################################################################
# AF_INPUTS - functions for modifying inputs before passing to alphafold
############################################################################
class _af_inputs:
  def _get_seq(self, params, opt, aux, key):
    '''get sequence features'''
    seq = {"input":params["seq"]}
    seq_shape = params["seq"].shape

    # shuffle msa (randomly pick which sequence is query)
    if self._num > 1:
      n = jax.random.randint(key,[],0,seq_shape[0])
      seq["input"] = seq["input"].at[0].set(seq["input"][n]).at[n].set(seq["input"][0])

    # straight-through/reparameterization
    seq["logits"] = seq["input"] * opt["alpha"] + opt["bias"]
    seq["soft"] = jax.nn.softmax(seq["logits"] / opt["temp"])
    seq["hard"] = jax.nn.one_hot(seq["soft"].argmax(-1), 20)
    seq["hard"] = jax.lax.stop_gradient(seq["hard"] - seq["soft"]) + seq["soft"]

    # create pseudo sequence
    seq["pseudo"] = opt["soft"] * seq["soft"] + (1-opt["soft"]) * seq["input"]
    seq["pseudo"] = opt["hard"] * seq["hard"] + (1-opt["hard"]) * seq["pseudo"]
    
    # fix sequence for partial hallucination
    if self.protocol == "partial":
      seq_ref = jax.nn.one_hot(self._wt_aatype,20)
      fix_seq = lambda x:jnp.where(opt["fix_seq"],x.at[:,opt["pos"],:].set(seq_ref),x)
      seq = jax.tree_map(fix_seq,seq)

    aux.update({"seq":seq, "seq_pseudo":seq["pseudo"]})
    
    # protocol specific modifications to seq features
    if self.protocol == "binder":
      # concatenate target and binder sequence
      seq_target = jax.nn.one_hot(self._batch["aatype"][:self._target_len],20)
      seq_target = jnp.broadcast_to(seq_target,(self._num, *seq_target.shape))
      seq = jax.tree_map(lambda x:jnp.concatenate([seq_target,x],1), seq)
      
    if self.protocol in ["fixbb","hallucination"] and self._copies > 1:
      seq = jax.tree_map(lambda x:expand_copies(x, self._copies, self._args["block_diag"]), seq)

    return seq

  def _update_template(self, inputs, opt, key):
    ''''dynamically update template features'''
    if self.protocol in ["partial","fixbb","binder"]:      
      L = self._batch["aatype"].shape[0]
      
      if self.protocol in ["partial","fixbb"]:
        aatype = jnp.zeros(L)
        template_aatype = jnp.broadcast_to(opt["template"]["aatype"],(L,))
      
      if self.protocol == "binder":
        if self._redesign:
          aatype = jnp.asarray(self._batch["aatype"])
          aatype = aatype.at[self._target_len:].set(0)
          template_aatype = aatype.at[self._target_len:].set(opt["template"]["aatype"])
        else:
          aatype = template_aatype = self._batch["aatype"]
      
      # get pseudo-carbon-beta coordinates (carbon-alpha for glycine)
      pb, pb_mask = model.modules.pseudo_beta_fn(aatype,
                                                 self._batch["all_atom_positions"],
                                                 self._batch["all_atom_mask"])
      
      # define template features
      template_feats = {"template_aatype": template_aatype,
                        "template_all_atom_positions": self._batch["all_atom_positions"],
                        "template_all_atom_masks": self._batch["all_atom_mask"],
                        "template_pseudo_beta": pb,
                        "template_pseudo_beta_mask": pb_mask}
      
      # protocol specific template injection
      if self.protocol == "fixbb":
        for k,v in template_feats.items():
          inputs[k] = inputs[k].at[:,0].set(v)            
          if k == "template_all_atom_masks":
            inputs[k] = inputs[k].at[:,0,:,5:].set(0)
      
      if self.protocol == "binder":
        n = self._target_len
        for k,v in template_feats.items():
          inputs[k] = inputs[k].at[:,0,:n].set(v[:n])            
          inputs[k] = inputs[k].at[:,-1,n:].set(v[n:])
          if k == "template_all_atom_masks":
            inputs[k] = inputs[k].at[:,-1,n:,5:].set(0)       
      
      if self.protocol == "partial":
        p = opt["pos"]
        for k,v in template_feats.items():
          inputs[k] = inputs[k].at[:,0,p].set(v)
          # if k == "template_aatype": v = jax.nn.one_hot(v,22)
          # inputs[k] = jnp.einsum("ij,i...->j...",p,v)[None,None]
          if k == "template_all_atom_masks":
            inputs[k] = inputs[k].at[:,0,:,5:].set(0)

    # dropout template input features
    L = inputs["template_aatype"].shape[2]
    n = self._target_len if self.protocol == "binder" else 0
    pos_mask = jax.random.bernoulli(key, 1-opt["template"]["dropout"],(L,))
    inputs["template_all_atom_masks"] = inputs["template_all_atom_masks"].at[:,:,n:].multiply(pos_mask[n:,None])
    inputs["template_pseudo_beta_mask"] = inputs["template_pseudo_beta_mask"].at[:,:,n:].multiply(pos_mask[n:])

def update_seq(seq, inputs, seq_1hot=None, seq_pssm=None, msa_input=None):
  '''update the sequence features'''
  
  if seq_1hot is None: seq_1hot = seq 
  if seq_pssm is None: seq_pssm = seq
  
  seq_1hot = jnp.pad(seq_1hot,[[0,0],[0,0],[0,22-seq_1hot.shape[-1]]])
  seq_pssm = jnp.pad(seq_pssm,[[0,0],[0,0],[0,22-seq_pssm.shape[-1]]])
  
  msa_feat = jnp.zeros_like(inputs["msa_feat"]).at[...,0:22].set(seq_1hot).at[...,25:47].set(seq_pssm)
  if seq.ndim == 3:
    target_feat = jnp.zeros_like(inputs["target_feat"]).at[...,1:21].set(seq[0,...,:20])
  else:
    target_feat = jnp.zeros_like(inputs["target_feat"]).at[...,1:21].set(seq[...,:20])
    
  inputs.update({"target_feat":target_feat,"msa_feat":msa_feat})

def update_aatype(aatype, inputs):
  if jnp.issubdtype(aatype.dtype, jnp.integer):
    inputs.update({"aatype":aatype,
                   "atom14_atom_exists":residue_constants.restype_atom14_mask[aatype],
                   "atom37_atom_exists":residue_constants.restype_atom37_mask[aatype],
                   "residx_atom14_to_atom37":residue_constants.restype_atom14_to_atom37[aatype],
                   "residx_atom37_to_atom14":residue_constants.restype_atom37_to_atom14[aatype]})
  else:
    restype_atom14_to_atom37 = jax.nn.one_hot(residue_constants.restype_atom14_to_atom37,37)
    restype_atom37_to_atom14 = jax.nn.one_hot(residue_constants.restype_atom37_to_atom14,14)
    inputs.update({"aatype":aatype,
                   "atom14_atom_exists":jnp.einsum("...a,am->...m", aatype, residue_constants.restype_atom14_mask),
                   "atom37_atom_exists":jnp.einsum("...a,am->...m", aatype, residue_constants.restype_atom37_mask),
                   "residx_atom14_to_atom37":jnp.einsum("...a,abc->...bc", aatype, restype_atom14_to_atom37),
                   "residx_atom37_to_atom14":jnp.einsum("...a,abc->...bc", aatype, restype_atom37_to_atom14)})

def expand_copies(x, copies, block_diag=True):
  '''
  given msa (N,L,20) expand to (N*(1+copies),L*copies,22) if block_diag else (N,L*copies,22)
  '''
  if x.shape[-1] < 22:
    x = jnp.pad(x,[[0,0],[0,0],[0,22-x.shape[-1]]])
  x = jnp.tile(x,[1,copies,1])
  if copies > 1 and block_diag:
    L = x.shape[1]
    sub_L = L // copies
    y = x.reshape((-1,1,copies,sub_L,22))
    block_diag_mask = jnp.expand_dims(jnp.eye(copies),(0,3,4))
    seq = block_diag_mask * y
    gap_seq = (1-block_diag_mask) * jax.nn.one_hot(jnp.repeat(21,sub_L),22)  
    y = (seq + gap_seq).swapaxes(0,1).reshape(-1,L,22)
    return jnp.concatenate([x,y],0)
  else:
    return x