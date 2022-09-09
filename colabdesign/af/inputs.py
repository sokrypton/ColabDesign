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
  def _get_seq(self, inputs, params, opt, aux, key):
    '''get sequence features'''
    seq = soft_seq(params["seq"], opt, key)
    if "fix_pos" in opt:
      if "pos" in self.opt:
        seq_ref = jax.nn.one_hot(self._wt_aatype_sub,20)
        p = opt["pos"][opt["fix_pos"]]
        fix_seq = lambda x:x.at[...,p,:].set(seq_ref)
      else:
        seq_ref = jax.nn.one_hot(self._wt_aatype,20)
        p = opt["fix_pos"]
        fix_seq = lambda x:x.at[...,p,:].set(seq_ref[...,p,:])
      seq = jax.tree_map(fix_seq, seq)
    aux.update({"seq":seq, "seq_pseudo":seq["pseudo"]})
    
    # protocol specific modifications to seq features
    if self.protocol == "binder":
      # concatenate target and binder sequence
      seq_target = jax.nn.one_hot(inputs["batch"]["aatype"][:self._target_len],20)
      seq_target = jnp.broadcast_to(seq_target,(self._num, *seq_target.shape))
      seq = jax.tree_map(lambda x:jnp.concatenate([seq_target,x],1), seq)
      
    if self.protocol in ["fixbb","hallucination","partial"] and self._args["copies"] > 1:
      seq = jax.tree_map(lambda x:expand_copies(x, self._args["copies"], self._args["block_diag"]), seq)

    return seq

  def _update_template(self, inputs, opt, key):
    ''''dynamically update template features'''
    
    o = opt["template"]

    # enable templates
    inputs["template_mask"] = inputs["template_mask"].at[:].set(1)

    batch = inputs["batch"]
    if self.protocol in ["partial","fixbb","binder"]:

      L = batch["aatype"].shape[0]
      
      # decide which position to remove sequence and/or sidechains
      rm = jnp.logical_or(o["rm_seq"],o["rm_sc"])
      rm_seq = jnp.full(L,o["rm_seq"])
      rm_sc  = jnp.full(L,rm)
      if self.protocol == "binder":
        rm_seq = rm_seq.at[:self._target_len].set(False)
        rm_sc  = rm_sc.at[:self._target_len].set(False)

      # aatype = is used to define template's CB coordinates (CA in case of glycine)
      # template_aatype = is used as template's sequence
      aatype = jnp.where(rm_seq,0,batch["aatype"])
      template_aatype = jnp.where(rm_seq,21,batch["aatype"])
                    
      # get pseudo-carbon-beta coordinates (carbon-alpha for glycine)
      cb, cb_mask = model.modules.pseudo_beta_fn(aatype, batch["all_atom_positions"], batch["all_atom_mask"])
      
      # define template features
      template_feats = {"template_aatype": template_aatype,
                        "template_all_atom_positions": batch["all_atom_positions"],
                        "template_all_atom_mask": batch["all_atom_mask"],
                        "template_pseudo_beta": cb,
                        "template_pseudo_beta_mask": cb_mask}

      # inject template features
      if self.protocol == "partial":
        pos = opt["pos"]
        if self._args["repeat"] or self._args["homooligomer"]:
          C,L = self._args["copies"], self._len
          pos = (jnp.repeat(pos,C).reshape(-1,C) + jnp.arange(C) * L).T.flatten()

      for k,v in template_feats.items():
        if self.protocol == "partial":
          inputs[k] = inputs[k].at[0,pos].set(v)
        else:
          inputs[k] = inputs[k].at[0].set(v)
        
        # remove sidechains (mask anything beyond CB)
        if k == "template_all_atom_mask":
          if self.protocol == "partial":
            inputs[k] = inputs[k].at[:,pos,5:].set(jnp.where(rm_sc[:,None],0,inputs[k][:,pos,5:]))
          else:
            inputs[k] = inputs[k].at[:,:,5:].set(jnp.where(rm_sc[:,None],0,inputs[k][:,:,5:]))

    # dropout template input features
    L = inputs["template_aatype"].shape[1]
    n = self._target_len if self.protocol == "binder" else 0
    pos_mask = jax.random.bernoulli(key, 1-o["dropout"],(L,))
    inputs["template_all_atom_mask"] = inputs["template_all_atom_mask"].at[:,n:].multiply(pos_mask[n:,None])
    inputs["template_pseudo_beta_mask"] = inputs["template_pseudo_beta_mask"].at[:,n:].multiply(pos_mask[n:])

def update_seq(seq, inputs, seq_1hot=None, seq_pssm=None, mlm=None):
  '''update the sequence features'''
  
  if seq_1hot is None: seq_1hot = seq 
  if seq_pssm is None: seq_pssm = seq
  
  seq_1hot = jnp.pad(seq_1hot,[[0,0],[0,0],[0,22-seq_1hot.shape[-1]]])
  seq_pssm = jnp.pad(seq_pssm,[[0,0],[0,0],[0,22-seq_pssm.shape[-1]]])
  
  msa_feat = jnp.zeros_like(inputs["msa_feat"]).at[...,0:22].set(seq_1hot).at[...,25:47].set(seq_pssm)

  if mlm is not None:    
    X = jax.nn.one_hot(22,23)
    X = jnp.zeros(msa_feat.shape[-1]).at[...,:23].set(X).at[...,25:48].set(X)
    msa_feat = jnp.where(mlm[None,:,None],X,msa_feat)
    
  inputs.update({"msa_feat":msa_feat})

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
  given msa (N,L,20) expand to (1+N*copies,L*copies,22) if block_diag else (N,L*copies,22)
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
    return jnp.concatenate([x[:1],y],0)
  else:
    return x