import tensorflow as tf
tf.config.set_visible_devices([], 'GPU')

import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.shared.protein import *
from colabdesign.af.alphafold.common import protein, residue_constants
from colabdesign.af.alphafold.model import model, folding, all_atom

####################
# confidence metrics
####################
def get_plddt(outputs):
  logits = outputs["predicted_lddt"]["logits"]
  num_bins = logits.shape[-1]
  bin_width = 1.0 / num_bins
  bin_centers = jnp.arange(start=0.5 * bin_width, stop=1.0, step=bin_width)
  probs = jax.nn.softmax(logits, axis=-1)
  return jnp.sum(probs * bin_centers[None, :], axis=-1)

def get_pae(outputs):
  prob = jax.nn.softmax(outputs["predicted_aligned_error"]["logits"],-1)
  breaks = outputs["predicted_aligned_error"]["breaks"]
  step = breaks[1]-breaks[0]
  bin_centers = breaks + step/2
  bin_centers = jnp.append(bin_centers,bin_centers[-1]+step)
  return (prob*bin_centers).sum(-1)

####################
# loss functions
####################
def get_rmsd_loss_w(batch, outputs, copies=1):
  weights = batch["all_atom_mask"][:,1]
  true = batch["all_atom_positions"][:,1,:]
  pred = outputs["structure_module"]["final_atom_positions"][:,1,:]
  if copies == 1:
    return jnp_rmsd_w(true, pred, weights)
  else:
    # TODO add support for weights
    I = copies - 1
    L = true.shape[0] // copies
    p = true - true[:L].mean(0)
    q = pred - pred[:L].mean(0)
    p = p @ _np_kabsch(p[:L], q[:L])
    rm = jnp.square(p[:L]-q[:L]).sum(-1).mean()
    p,q = p[L:].reshape(I,1,L,-1),q[L:].reshape(1,I,L,-1)
    rm += jnp.square(p-q).sum(-1).mean(-1).min(-1).sum()
    return jnp.sqrt(rm / copies)

def get_rmsd_loss(batch, outputs):
  true = batch["all_atom_positions"][:,1,:]
  pred = outputs["structure_module"]["final_atom_positions"][:,1,:]
  return _np_rmsd(true,pred)

def _distogram_log_loss(logits, bin_edges, batch, num_bins, copies=1, return_errors=False):
  """Log loss of a distogram."""
  pos,mask = batch['pseudo_beta'],batch['pseudo_beta_mask']
  sq_breaks = jnp.square(bin_edges)
  dist2 = jnp.square(pos[:,None] - pos[None,:]).sum(-1,keepdims=True)
  true_bins = jnp.sum(dist2 > sq_breaks, axis=-1)
  true = jax.nn.one_hot(true_bins, num_bins)

  if copies == 1:
    errors = -(true * jax.nn.log_softmax(logits)).sum(-1)    
    sq_mask = mask[:,None] * mask[None,:]
    avg_error = (errors * sq_mask).sum()/(1e-6 + sq_mask.sum())
    if return_errors:
      return {"errors": errors, "loss": avg_error}
    else:
      return avg_error
  else:
    # TODO add support for masks
    L = pos.shape[0] // copies
    I = copies - 1
    true_, pred_ = true[:L,:L], logits[:L,:L]
    intra_errors = -(true_ * jax.nn.log_softmax(pred_)).sum(-1)
    avg_error = intra_errors.mean()

    true_, pred_ = true[:L,L:], logits[:L,L:]
    true_, pred_ = true_.reshape(L,I,1,L,-1), pred_.reshape(L,1,I,L,-1)
    inter_errors = -(true_ * jax.nn.log_softmax(pred_)).sum(-1)
    avg_error += inter_errors.mean((0,-1)).min(-1).sum()
    
    if return_errors:
      return {"intra_errors": intra_errors,
              "inter_errors": inter_errors,
              "loss": avg_error/copies}
    else:
      return avg_error/copies

def get_dgram_loss(batch, outputs, model_config, logits=None, aatype=None, copies=1, return_errors=False):
  # get cb features (ca in case of glycine)
  if aatype is None: aatype = batch["aatype"]
  pb, pb_mask = model.modules.pseudo_beta_fn(aatype,
                                             batch["all_atom_positions"],
                                             batch["all_atom_mask"])
  if logits is None: logits = outputs["distogram"]["logits"]
  dgram_loss = _distogram_log_loss(logits,
                                   outputs["distogram"]["bin_edges"],
                                   batch={"pseudo_beta":pb,"pseudo_beta_mask":pb_mask},
                                   num_bins=model_config.model.heads.distogram.num_bins,
                                   copies=copies, return_errors=return_errors)
  return dgram_loss

def get_fape_loss(batch, outputs, model_config, use_clamped_fape=False):
  sub_batch = jax.tree_map(lambda x: x, batch)
  sub_batch["use_clamped_fape"] = use_clamped_fape
  loss = {"loss":0.0}    
  folding.backbone_loss(loss, sub_batch, outputs["structure_module"], model_config.model.heads.structure_module)
  return loss["loss"]

def get_sc_rmsd(true_pos, pred_pos, aa_ident, atoms_to_exclude=None):
  if atoms_to_exclude is None: atoms_to_exclude = ["N","C","O"]
  
  # collect atom indices
  idx,idx_alt = [],[]
  for n,a in enumerate(aa_ident):
    aa = idx_to_resname[a]
    atoms = set(residue_constants.residue_atoms[aa])
    atoms14 = residue_constants.restype_name_to_atom14_names[aa]
    swaps = residue_constants.residue_atom_renaming_swaps.get(aa,{})
    swaps.update({v:k for k,v in swaps.items()})
    for atom in atoms.difference(atoms_to_exclude):
      idx.append(n * 14 + atoms14.index(atom))
      if atom in swaps:
        idx_alt.append(n * 14 + atoms14.index(swaps[atom]))
      else:
        idx_alt.append(idx[-1])
  idx, idx_alt = np.asarray(idx), np.asarray(idx_alt)

  # select atoms
  T, P = true_pos.reshape(-1,3)[idx], pred_pos.reshape(-1,3)[idx]

  # select non-ambigious atoms
  non_amb = idx == idx_alt
  t, p = T[non_amb], P[non_amb]

  # align non-ambigious atoms
  aln = _np_kabsch(t-t.mean(0), p-p.mean(0))
  T,P = (T-t.mean(0)) @ aln, P-p.mean(0)
  P_alt = pred_pos.reshape(-1,3)[idx_alt]-p.mean(0)

  # compute rmsd
  msd = jnp.minimum(jnp.square(T-P).sum(-1),jnp.square(T-P_alt).sum(-1)).mean()
  return jnp.sqrt(msd + 1e-8)

def get_6D_loss(batch, outputs, **kwargs):
  true = batch["all_atom_positions"]
  pred = outputs["structure_module"]["final_atom_positions"]
  mask = batch["all_atom_mask"]
  return _np_get_6D_loss(true, pred, mask, **kwargs)

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
  