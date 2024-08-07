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
    pssm = jnp.where(opt["pssm_hard"], seq["hard"], seq["pseudo"])
    if self._args["use_mlm"]:
      shape = seq["pseudo"].shape[:2]
      aux["mlm_mask"] = jax.random.bernoulli(k2,opt["mlm_dropout"],shape)
      update_seq(seq, inputs, seq_pssm=pssm, mlm=aux["mlm_mask"], mask_target=self._args["mask_target"])
    else:
      update_seq(seq, inputs, seq_pssm=pssm)
    
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

  def _update_template(self, inputs, key):
    ''''dynamically update template features''' 
    if "batch" in inputs:
      batch, opt = inputs["batch"], inputs["opt"]

      # enable templates
      inputs["template_mask"] = inputs["template_mask"].at[0].set(1)
      L = batch["aatype"].shape[0]
      
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
        # get pseudo-carbon-beta coordinates (carbon-alpha for glycine)
        # aatype = is used to define template's CB coordinates (CA in case of glycine)
        cb, cb_mask = model.modules.pseudo_beta_fn(
          jnp.where(rm_seq,0,batch["aatype"]),
          batch["all_atom_positions"],
          batch["all_atom_mask"])
        template_feats.update({"template_pseudo_beta":        cb,
                               "template_pseudo_beta_mask":   cb_mask,
                               "template_all_atom_positions": batch["all_atom_positions"],
                               "template_all_atom_mask":      batch["all_atom_mask"]})

      # inject template features
      if self.protocol == "partial":
        pos = opt["pos"]
        if self._args["repeat"] or self._args["homooligomer"]:
          C,L = self._args["copies"], self._len
          pos = (jnp.repeat(pos,C).reshape(-1,C) + jnp.arange(C) * L).T.flatten()

      for k,v in template_feats.items():
        if self.protocol == "partial":
          if k in ["template_dgram"]:
            inputs[k] = inputs[k].at[0,pos[:,None],pos[None,:]].set(v)
          else:
            inputs[k] = inputs[k].at[0,pos].set(v)
        else:
          inputs[k] = inputs[k].at[0].set(v)
        
        # remove sidechains (mask anything beyond CB)
        if k in ["template_all_atom_mask"]:
          if self.protocol == "partial":
            inputs[k] = inputs[k].at[:,pos,5:].set(jnp.where(rm_sc[:,None],0,inputs[k][:,pos,5:]))
            inputs[k] = inputs[k].at[:,pos].set(jnp.where(rm[:,None],0,inputs[k][:,pos]))
          else:
            inputs[k] = inputs[k].at[...,5:].set(jnp.where(rm_sc[:,None],0,inputs[k][...,5:]))
            inputs[k] = jnp.where(rm[:,None],0,inputs[k])

def update_seq(seq, inputs, seq_1hot=None, seq_pssm=None, mlm=None, mask_target=False, use_jax=True):
  '''update the sequence features'''  

  _np = jnp if use_jax else np
  if isinstance(seq, dict):
    if seq_1hot is None: seq_1hot = seq["pseudo"] 
    if seq_pssm is None: seq_pssm = seq["pseudo"]
  else:
    if seq_1hot is None: seq_1hot = seq
    if seq_pssm is None: seq_pssm = seq
    
  target_feat = seq_1hot[0,:,:20]

  seq_1hot = _np.pad(seq_1hot,[[0,0],[0,0],[0,22-seq_1hot.shape[-1]]])
  seq_pssm = _np.pad(seq_pssm,[[0,0],[0,0],[0,22-seq_pssm.shape[-1]]])
  if use_jax:
    msa_feat = jnp.zeros_like(inputs["msa_feat"]).at[...,0:22].set(seq_1hot).at[...,25:47].set(seq_pssm)
  else:
    msa_feat = np.zeros_like(inputs["msa_feat"])
    msa_feat[...,0:22] = seq_1hot
    msa_feat[...,25:47] = seq_pssm


  # masked language modeling (randomly mask positions)
  if mlm is not None:
    X = _np.eye(23)[22]
    if use_jax:
      Y = jnp.zeros(msa_feat.shape[-1]).at[...,:23].set(X).at[...,25:48].set(X)
    else:
      Y = np.zeros(msa_feat.shape[-1])
      Y[...,:23] = X 
      Y[...,25:48] = X

    msa_feat = _np.where(mlm[...,None],Y,msa_feat)
    seq["pseudo"] = _np.where(mlm[...,None],X[:seq["pseudo"].shape[-1]],seq["pseudo"])
    if mask_target:
      target_feat = _np.where(mlm[0,:,None],0,target_feat)
    
  inputs.update({"msa_feat":msa_feat, "target_feat":target_feat})

def update_aatype(seq, inputs, use_jax=True):
  _np = jnp if use_jax else np
  if isinstance(seq,dict):
    aatype = seq["pseudo"][0].argmax(-1)
  else:
    aatype = seq[0].argmax(-1)
  r = residue_constants
  a = {"atom14_atom_exists":r.restype_atom14_mask,
       "atom37_atom_exists":r.restype_atom37_mask,
       "residx_atom14_to_atom37":r.restype_atom14_to_atom37,
       "residx_atom37_to_atom14":r.restype_atom37_to_atom14}
  mask = inputs["seq_mask"][:,None]
  inputs.update(jax.tree_map(lambda x:_np.where(mask,_np.asarray(x)[aatype],0),a))
  inputs["aatype"] = aatype

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