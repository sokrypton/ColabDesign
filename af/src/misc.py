import jax
import jax.numpy as jnp
import tensorflow as tf
tf.config.set_visible_devices([], 'GPU')

import numpy as np
from alphafold.common import protein, residue_constants
from alphafold.model import model, folding, all_atom

#########################
# rmsd
#########################
def jnp_rmsdist(true, pred):
  return _np_rmsdist(true, pred)

def jnp_rmsd(true, pred, add_dist=False):
  rmsd = _np_rmsd(true, pred)
  if add_dist: rmsd = (rmsd + _np_rmsdist(true, pred))/2
  return rmsd

def jnp_kabsch_w(a, b, weights):
  return _np_kabsch(a * weights[:,None], b)

def jnp_rmsd_w(true, pred, weights):
  p = true - (true * weights[:,None]).sum(0,keepdims=True)/weights.sum()
  q = pred - (pred * weights[:,None]).sum(0,keepdims=True)/weights.sum()
  p = p @ _np_kabsch(p * weights[:,None], q)
  return jnp.sqrt((weights*jnp.square(p-q).sum(-1)).sum()/weights.sum() + 1e-8)

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
def get_rmsd_loss(batch, outputs):
  true = batch["all_atom_positions"][:,1,:]
  pred = outputs["structure_module"]["final_atom_positions"][:,1,:]
  return _np_rmsd(true,pred)

def _distogram_log_loss(logits, bin_edges, batch, num_bins, copies=1):
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
    return avg_error
  else:
    # TODO add support for masks
    L = pos.shape[0] // copies
    I = copies - 1
    true_, pred_ = true[:L,:L], logits[:L,:L]
    errors = -(true_ * jax.nn.log_softmax(pred_)).sum(-1)
    avg_error = errors.mean()

    true_, pred_ = true[:L,L:], logits[:L,L:]
    true_, pred_ = true_.reshape(L,I,1,L,-1), pred_.reshape(L,1,I,L,-1)
    errors = -(true_ * jax.nn.log_softmax(pred_)).sum(-1)
    avg_error += errors.mean((0,-1)).min(-1).sum()

    return avg_error / copies

def get_dgram_loss(batch, outputs, model_config, logits=None, copies=1):
  # get cb features (ca in case of glycine)
  pb, pb_mask = model.modules.pseudo_beta_fn(batch["aatype"],
                                             batch["all_atom_positions"],
                                             batch["all_atom_mask"])
  if logits is None: logits = outputs["distogram"]["logits"]
  dgram_loss = _distogram_log_loss(logits,
                                   outputs["distogram"]["bin_edges"],
                                   batch={"pseudo_beta":pb,"pseudo_beta_mask":pb_mask},
                                   num_bins=model_config.model.heads.distogram.num_bins,
                                   copies=copies)
  return dgram_loss

def get_fape_loss(batch, outputs, model_config, use_clamped_fape=False):
  sub_batch = jax.tree_map(lambda x: x, batch)
  sub_batch["use_clamped_fape"] = use_clamped_fape
  loss = {"loss":0.0}    
  folding.backbone_loss(loss, sub_batch, outputs["structure_module"], model_config.model.heads.structure_module)
  return loss["loss"]

####################
# loss functions (restricted to idx and/or sidechains)
####################
def get_dgram_loss_idx(batch, outputs, idx, model_config):
  idx_ref = batch["idx"]
  pb, pb_mask = model.modules.pseudo_beta_fn(batch["aatype"][idx_ref],
                                             batch["all_atom_positions"][idx_ref],
                                             batch["all_atom_mask"][idx_ref])
  
  dgram_loss = model.modules._distogram_log_loss(outputs["distogram"]["logits"][:,idx][idx,:],
                                                 outputs["distogram"]["bin_edges"],
                                                 batch={"pseudo_beta":pb,"pseudo_beta_mask":pb_mask},
                                                 num_bins=model_config.model.heads.distogram.num_bins)
  return dgram_loss["loss"]

def get_fape_loss_idx(batch, outputs, idx, model_config, backbone=False, sidechain=True, use_clamped_fape=False):
  idx_ref = batch["idx"]
  
  sub_batch = batch.copy()
  sub_batch.pop("idx")
  sub_batch = jax.tree_map(lambda x: x[idx_ref,...],sub_batch)
  sub_batch["use_clamped_fape"] = use_clamped_fape
  
  value = jax.tree_map(lambda x: x, outputs["structure_module"])
  loss = {"loss":0.0}
  
  if sidechain:
    value.update(folding.compute_renamed_ground_truth(sub_batch, value['final_atom14_positions'][idx,...]))
    value['sidechains']['frames'] = jax.tree_map(lambda x: x[:,idx,:], value["sidechains"]["frames"])
    value['sidechains']['atom_pos'] = jax.tree_map(lambda x: x[:,idx,:], value["sidechains"]["atom_pos"])
    loss.update(folding.sidechain_loss(sub_batch, value, model_config.model.heads.structure_module))
  
  if backbone:
    value["traj"] = value["traj"][...,idx,:]
    folding.backbone_loss(loss, sub_batch, value, model_config.model.heads.structure_module)

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

def get_sidechain_rmsd_idx(batch, outputs, idx, model_config, include_ca=True):
  idx_ref = batch["idx"]
  true_aa_idx = batch["aatype"][idx_ref]
  true_pos = all_atom.atom37_to_atom14(batch["all_atom_positions"],batch)[idx_ref,:,:]
  pred_pos = outputs["structure_module"]["final_atom14_positions"][idx,:,:]
  bb_atoms_to_exclude = ["N","C","O"] if include_ca else ["N","CA","C","O"]
  
  return get_sc_rmsd(true_pos, pred_pos, true_aa_idx, bb_atoms_to_exclude)

#################################################################################
#################################################################################
#################################################################################

def _np_len_pw(x, use_jax=True):
  '''compute pairwise distance'''
  _np = jnp if use_jax else np

  x_norm = _np.square(x).sum(-1)
  xx = _np.einsum("...ia,...ja->...ij",x,x)
  sq_dist = x_norm[...,:,None] + x_norm[...,None,:] - 2 * xx

  # due to precision errors the values can sometimes be negative
  if use_jax: sq_dist = jax.nn.relu(sq_dist)
  else: sq_dist[sq_dist < 0] = 0

  # return euclidean pairwise distance matrix
  return _np.sqrt(sq_dist + 1e-8)

def _np_rmsdist(true, pred, use_jax=True):
  '''compute RMSD of distance matrices'''
  _np = jnp if use_jax else np
  t = _np_len_pw(true, use_jax=use_jax)
  p = _np_len_pw(pred, use_jax=use_jax)
  return _np.sqrt(_np.square(t-p).mean() + 1e-8)

def _np_kabsch(a, b, return_v=False, use_jax=True):
  '''get alignment matrix for two sets of coodinates'''
  _np = jnp if use_jax else np
  ab = a.swapaxes(-1,-2) @ b
  u, s, vh = _np.linalg.svd(ab, full_matrices=False)
  flip = _np.linalg.det(u @ vh) < 0
  u_ = _np.where(flip, -u[...,-1].T, u[...,-1].T).T
  if use_jax: u = u.at[...,-1].set(u_)
  else: u[...,-1] = u_
  return u if return_v else (u @ vh)

def _np_rmsd(true, pred, use_jax=True):
  '''compute RMSD of coordinates after alignment'''
  _np = jnp if use_jax else np
  p = true - true.mean(-2,keepdims=True)
  q = pred - pred.mean(-2,keepdims=True)
  p = p @ _np_kabsch(p, q, use_jax=use_jax)
  return _np.sqrt(_np.square(p-q).sum(-1).mean(-1) + 1e-8)

def _np_norm(x, axis=-1, keepdims=True, eps=1e-8, use_jax=True):
  '''compute norm of vector'''
  _np = jnp if use_jax else np
  return _np.sqrt(_np.square(x).sum(axis,keepdims=keepdims) + 1e-8)
  
def _np_len(a, b, use_jax=True):
  '''given coordinates a-b, return length or distance'''
  return _np_norm(a-b, use_jax=use_jax)

def _np_ang(a, b, c, use_acos=False, use_jax=True):
  '''given coordinates a-b-c, return angle'''  
  _np = jnp if use_jax else np
  norm = lambda x: _np_norm(x, use_jax=use_jax)
  ba, bc = b-a, b-c
  cos_ang = (ba * bc).sum(-1,keepdims=True) / (norm(ba) * norm(bc))
  # note the derivative at acos(-1 or 1) is inf, to avoid nans we use cos(ang)
  if use_acos: return _np.arccos(cos_ang)
  else: return cos_ang
  
def _np_dih(a, b, c, d, use_atan2=False, standardize=False, use_jax=True):
  '''given coordinates a-b-c-d, return dihedral'''
  _np = jnp if use_jax else np
  normalize = lambda x: x/_np_norm(x, use_jax=use_jax)
  ab, bc, cd = normalize(a-b), normalize(b-c), normalize(c-d)
  n1,n2 = _np.cross(ab, bc), _np.cross(bc, cd)
  sin_ang = (_np.cross(n1, bc) * n2).sum(-1,keepdims=True)
  cos_ang = (n1 * n2).sum(-1,keepdims=True)
  if use_atan2:
    return _np.arctan2(sin_ang, cos_ang)
  else:
    angs = _np.concatenate([sin_ang, cos_ang],-1)
    if standardize: return normalize(angs)
    else: return angs

def _np_extend(a,b,c, L,A,D, use_jax=True):
  '''
  given coordinates a-b-c,
  c-d (L)ength, b-c-d (A)ngle, and a-b-c-d (D)ihedral
  return 4th coordinate d
  '''
  _np = jnp if use_jax else np
  normalize = lambda x: x/_np_norm(x, use_jax=use_jax)
  bc = normalize(b-c)
  n = normalize(_np.cross(b-a, bc))
  return c + sum([L * _np.cos(A) * bc,
                  L * _np.sin(A) * _np.cos(D) * _np.cross(n, bc),
                  L * _np.sin(A) * _np.sin(D) * -n])

def _np_get_cb(N,CA,C, use_jax=True):
  '''compute CB placement from N, CA, C'''
  return _np_extend(C, N, CA, 1.522, 1.927, -2.143, use_jax=use_jax)
  
def _np_get_6D(all_atom_positions, all_atom_mask=None, use_jax=True):
  '''get 6D features (see TrRosetta paper)'''

  # get CB coordinate
  atom_idx = {k:residue_constants.atom_order[k] for k in ["N","CA","C"]}
  out = {k:all_atom_positions[...,i,:] for k,i in atom_idx.items()}
  out["CB"] = _np_get_cb(**out, use_jax=use_jax)
  
  if all_atom_mask is not None:
    idx = np.fromiter(atom_idx.values(),int)
    out["CB_mask"] = all_atom_mask[...,idx].prod(-1)

  # get pairwise features
  N,A,B = (out[k] for k in ["N","CA","CB"])
  j = {"use_jax":use_jax}
  out.update({"dist":  _np_len_pw(B,**j),
              "phi":   _np_ang(A[...,:,None,:],B[...,:,None,:],B[...,None,:,:],**j),
              "omega": _np_dih(A[...,:,None,:],B[...,:,None,:],B[...,None,:,:],A[...,None,:,:],**j),
              "theta": _np_dih(N[...,:,None,:],A[...,:,None,:],B[...,:,None,:],B[...,None,:,:],**j),
              })
  return out

####################
# 6D loss (see TrRosetta paper)
####################
def _np_get_6D_loss(true, pred, mask=None, use_theta=True, use_dist=False, use_jax=True):
  _np = jnp if use_jax else np

  f = {"T":_np_get_6D(true, mask, use_jax=use_jax),
       "P":_np_get_6D(pred, use_jax=use_jax)}

  for k in f: f[k]["dist"] /= 10.0

  keys = ["omega","phi"]
  if use_theta: keys.append("theta")
  if use_dist: keys.append("dist")
  sq_diff = sum([_np.square(f["T"][k]-f["P"][k]).sum(-1) for k in keys])

  mask = _np.ones(true.shape[0]) if mask is None else f["T"]["CB_mask"]
  mask = mask[:,None] * mask[None,:]
  loss = (sq_diff * mask).sum((-1,-2)) / mask.sum((-1,-2))

  return _np.sqrt(loss + 1e-8).mean()
  
def get_6D_loss(batch, outputs, **kwargs):
  true = batch["all_atom_positions"]
  pred = outputs["structure_module"]["final_atom_positions"]
  mask = batch["all_atom_mask"]
  return _np_get_6D_loss(true, pred, mask, **kwargs)

#################################################################################
#################################################################################
#################################################################################

####################
# update sequence
####################
def soft_seq(seq_logits, temp=1.0, hard=True):
  seq_soft = jax.nn.softmax(seq_logits / temp)
  if hard:
    seq_hard = jax.nn.one_hot(seq_soft.argmax(-1),20)
    return jax.lax.stop_gradient(seq_hard - seq_soft) + seq_soft
  else:
    return seq_soft

def update_seq(seq, inputs, seq_1hot=None, seq_pssm=None, msa_input=None):
  '''update the sequence features'''
  
  if seq_1hot is None: seq_1hot = seq 
  if seq_pssm is None: seq_pssm = seq
  msa_feat = jnp.zeros_like(inputs["msa_feat"]).at[...,0:20].set(seq_1hot).at[...,25:45].set(seq_pssm)
  if seq.ndim == 3:
    target_feat = jnp.zeros_like(inputs["target_feat"]).at[...,1:21].set(seq[0])
  else:
    target_feat = jnp.zeros_like(inputs["target_feat"]).at[...,1:21].set(seq)
    
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
  
####################
# utils
####################

def pdb_to_string(pdb_file):
  lines = []
  for line in open(pdb_file,"r"):
    if line[:6] == "HETATM" and line[17:20] == "MSE":
      line = "ATOM  "+line[6:17]+"MET"+line[20:]
    if line[:4] == "ATOM":
      lines.append(line)
  return "".join(lines)

def save_pdb(outs, filename="tmp.pdb"):
  seq = outs["seq"].argmax(-1)
  while seq.ndim > 1: seq = seq[0]
  b_factors = np.zeros_like(outs["outputs"]['final_atom_mask'])
  p = protein.Protein(
        aatype=seq,
        atom_positions=outs["outputs"]["final_atom_positions"],
        atom_mask=outs["outputs"]['final_atom_mask'],
        residue_index=jnp.arange(len(seq))+1,
        b_factors=b_factors)
  pdb_lines = protein.to_pdb(p)
  with open(filename, 'w') as f:
    f.write(pdb_lines)

order_restype = {v: k for k, v in residue_constants.restype_order.items()}
idx_to_resname = dict((v,k) for k,v in residue_constants.resname_to_idx.items())

###########################
# MISC
###########################
jalview_color_list = {"Clustal":           ["#80a0f0","#f01505","#00ff00","#c048c0","#f08080","#00ff00","#c048c0","#f09048","#15a4a4","#80a0f0","#80a0f0","#f01505","#80a0f0","#80a0f0","#ffff00","#00ff00","#00ff00","#80a0f0","#15a4a4","#80a0f0"],
                      "Zappo":             ["#ffafaf","#6464ff","#00ff00","#ff0000","#ffff00","#00ff00","#ff0000","#ff00ff","#6464ff","#ffafaf","#ffafaf","#6464ff","#ffafaf","#ffc800","#ff00ff","#00ff00","#00ff00","#ffc800","#ffc800","#ffafaf"],
                      "Taylor":            ["#ccff00","#0000ff","#cc00ff","#ff0000","#ffff00","#ff00cc","#ff0066","#ff9900","#0066ff","#66ff00","#33ff00","#6600ff","#00ff00","#00ff66","#ffcc00","#ff3300","#ff6600","#00ccff","#00ffcc","#99ff00"],
                      "Hydrophobicity":    ["#ad0052","#0000ff","#0c00f3","#0c00f3","#c2003d","#0c00f3","#0c00f3","#6a0095","#1500ea","#ff0000","#ea0015","#0000ff","#b0004f","#cb0034","#4600b9","#5e00a1","#61009e","#5b00a4","#4f00b0","#f60009","#0c00f3","#680097","#0c00f3"],
                      "Helix Propensity":  ["#e718e7","#6f906f","#1be41b","#778877","#23dc23","#926d92","#ff00ff","#00ff00","#758a75","#8a758a","#ae51ae","#a05fa0","#ef10ef","#986798","#00ff00","#36c936","#47b847","#8a758a","#21de21","#857a85","#49b649","#758a75","#c936c9"],
                      "Strand Propensity": ["#5858a7","#6b6b94","#64649b","#2121de","#9d9d62","#8c8c73","#0000ff","#4949b6","#60609f","#ecec13","#b2b24d","#4747b8","#82827d","#c2c23d","#2323dc","#4949b6","#9d9d62","#c0c03f","#d3d32c","#ffff00","#4343bc","#797986","#4747b8"],
                      "Turn Propensity":   ["#2cd3d3","#708f8f","#ff0000","#e81717","#a85757","#3fc0c0","#778888","#ff0000","#708f8f","#00ffff","#1ce3e3","#7e8181","#1ee1e1","#1ee1e1","#f60909","#e11e1e","#738c8c","#738c8c","#9d6262","#07f8f8","#f30c0c","#7c8383","#5ba4a4"],
                      "Buried Index":      ["#00a35c","#00fc03","#00eb14","#00eb14","#0000ff","#00f10e","#00f10e","#009d62","#00d52a","#0054ab","#007b84","#00ff00","#009768","#008778","#00e01f","#00d52a","#00db24","#00a857","#00e619","#005fa0","#00eb14","#00b649","#00f10e"]}
