import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.shared.model import soft_seq
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.af.alphafold.model import model

############################################################################
# AF_INPUTS - functions for modifying inputs before passing to alphafold
############################################################################
class _af_inputs:
  def _get_seq(self, params, opt, aux, key):
    '''get sequence features'''
    seq = soft_seq(params["seq"], opt, key)
    if "pos" in opt and "fix_seq" in opt:
      seq_ref = jax.nn.one_hot(self._wt_aatype,20)
      p = opt["pos"]
      if self.protocol == "partial":
        fix_seq = lambda x:jnp.where(opt["fix_seq"],x.at[...,p,:].set(seq_ref),x)
      else:
        fix_seq = lambda x:jnp.where(opt["fix_seq"],x.at[...,p,:].set(seq_ref[...,p,:]),x)
      seq = jax.tree_map(fix_seq,seq)
    aux.update({"seq":seq, "seq_pseudo":seq["pseudo"]})
    
    # protocol specific modifications to seq features
    if self.protocol == "binder":
      # concatenate target and binder sequence
      seq_target = jax.nn.one_hot(self._batch["aatype"][:self._target_len],20)
      seq_target = jnp.broadcast_to(seq_target,(self._num, *seq_target.shape))
      seq = jax.tree_map(lambda x:jnp.concatenate([seq_target,x],1), seq)
      
    if self.protocol in ["fixbb","hallucination"] and self._args["copies"] > 1:
      seq = jax.tree_map(lambda x:expand_copies(x, self._args["copies"], self._args["block_diag"]), seq)

    return seq

  def _update_template(self, inputs, opt, key):
    ''''dynamically update template features'''

    # aatype = is used to define template's CB coordinates (CA in case of glycine)
    # template_aatype = is used as template's sequence

    if self.protocol in ["partial","fixbb","binder"]:      
      L = self._batch["aatype"].shape[0]
      
      if self.protocol in ["partial","fixbb"]:
        if self._args["rm_template_seq"]:
          aatype = jnp.zeros(L) 
          template_aatype = jnp.broadcast_to(opt["template"]["aatype"],(L,))
        else:
          template_aatype = aatype = self._batch["aatype"]
      
      if self.protocol == "binder":
        if self._args["redesign"] and self._args["rm_template_seq"]:
          aatype          = jnp.asarray(self._batch["aatype"])
          aatype          = aatype.at[self._target_len:].set(0)
          template_aatype = aatype.at[self._target_len:].set(opt["template"]["aatype"])
        
        else:
          # target=(template_seq=True)        
          template_aatype = aatype = self._batch["aatype"]
      
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
      for k,v in template_feats.items():
        if self.protocol == "binder":
          n = self._target_len
          inputs[k] = inputs[k].at[:,0,:n].set(v[:n])
          inputs[k] = inputs[k].at[:,-1,n:].set(v[n:])
        
        if self.protocol == "fixbb":
          inputs[k] = inputs[k].at[:,0].set(v)

        if self.protocol == "partial":
          inputs[k] = inputs[k].at[:,0,opt["pos"]].set(v)
        
        if k == "template_all_atom_masks" and self._args["rm_template_seq"]:
          if self.protocol == "binder":
            inputs[k] = inputs[k].at[:,-1,n:,5:].set(0)
          else:
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

def crop_feat(feat, pos, cfg, add_batch=True):  
  def find(x,k):
    i = []
    for j,y in enumerate(x):
      if y == k: i.append(j)
    return i
  if feat is None:
    return None
  else:
    shapes = cfg.data.eval.feat
    NUM_RES = "num residues placeholder"
    idx = {k:find(v,NUM_RES) for k,v in shapes.items()}

    new_feat = {}
    for k,v in feat.items():
      v_ = v.copy()
      if k in idx:
        for i in idx[k]:
          v_ = jnp.take(v_, pos, i + add_batch)
      new_feat[k] = v_
    return new_feat