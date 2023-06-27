import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.shared.utils import Key, copy_dict
from colabdesign.shared.protein import jnp_rmsd_w, _np_kabsch, _np_rmsd, _np_get_6D_loss
from colabdesign.af.alphafold.model import model, folding, all_atom, geometry
from colabdesign.af.alphafold.common import confidence, residue_constants
from colabdesign.af.alphafold.model import folding_multimer, all_atom_multimer

####################################################
# AF_LOSS - setup loss function
####################################################

class _af_loss:
  # protocol specific loss functions

  def _get_loss(self, inputs, outputs, aux):
    opt = inputs["opt"]
    self._loss_supervised(inputs, outputs, aux)
    self._loss_unsupervised(inputs, outputs, aux)

  def _loss_supervised(self, inputs, outputs, aux):
    a,opt = self._args, inputs["opt"]

    if self.protocol == "binder":
      # TODO: remove protocol specific losses
      aln = get_rmsd_loss(inputs, outputs, copies=a["copies"], include_L=False, L=self._target_len)
    else:
      aln = get_rmsd_loss(inputs, outputs, copies=a["copies"])
    if a["realign"]:
      aux["atom_positions"] = aln["align"](aux["atom_positions"]) * aux["atom_mask"][...,None]
    
    aux["losses"].update({
      "fape":      get_fape_loss(inputs, outputs, copies=a["copies"], clamp=opt["fape_cutoff"]),
      "dgram_cce": get_dgram_loss(inputs, outputs, copies=a["copies"], aatype=inputs["aatype"]),
      "rmsd":      aln["rmsd"],
    })
    if a["use_sidechains"]:
      if "fix_pos" in inputs:
        aux["losses"].update(get_sc_loss(inputs, outputs, cfg=self._cfg, mask_1d=inputs["fix_pos"]))    
      # TODO: add support for optimize_seq=False

  def _loss_unsupervised(self, inputs, outputs, aux, mask_1d=None, mask_2d=None):

    opt = inputs["opt"]

    # define masks
    seq_mask_1d = inputs["seq_mask"]
    seq_mask_2d = seq_mask_1d[:,None] * seq_mask_1d[None,:]
    intra_mask_2d = inputs["asym_id"][:,None] == inputs["asym_id"][None,:]
    m = lambda a,b,c: jnp.where(a,b,c)

    # partial masks
    mask_1d = m(seq_mask_1d,inputs["unsupervised"],0)
    mask_2d = m(seq_mask_2d,inputs["unsupervised_2d"],0)      
    masks_partial = {"mask_1d":mask_1d,              "mask_2d":mask_2d,
     "intra_masks": {"mask_1d":mask_2d_intra.max(1), "mask_2d":m(intra_mask_2d,mask_2d,0)},
     "inter_masks": {"mask_1d":mask_2d_inter.max(1), "mask_2d":m(intra_mask_2d,0,mask_2d)}}
    
    # full masks
    masks_full =    {"mask_1d":seq_mask_1d, "mask_2d":seq_mask_2d,
      "intra_masks":{"mask_1d":seq_mask_1d, "mask_2d":m(intra_mask_2d,seq_mask_2d,0)},
      "inter_masks":{"mask_1d":seq_mask_1d, "mask_2d":m(intra_mask_2d,0,seq_mask_2d)}}
    
    masks = jax.tree_map(lambda x,y: m(opt["partial_loss"],x,y), partial_masks, full_masks)
    if "hotspot" in inputs: masks["inter_masks"]["mask_1d"] = m(inputs["hotspot"],seq_mask_1d,0)
    aux["masks"] = masks

    # define losses
    losses = {
      "exp_res": get_exp_res_loss(outputs, mask_1d=masks["mask_1d"]),
      "plddt":   get_plddt_loss(outputs, mask_1d=masks["mask_1d"]),
      "pae":     get_pae_loss(outputs, mask_2d=masks["intra_masks"]["mask_2d"]),
      "con":     get_con_loss(inputs, outputs, opt["con"], **masks["intra_masks"]),
      "helix":   get_helix_loss(inputs, outputs, mask_2d=masks["mask_2d"]),
      "i_pae":   get_pae_loss(outputs, mask_2d=masks["inter_masks"]["mask_2d"]),
      "i_con":   get_con_loss(inputs, outputs, opt["i_con"], **masks["inter_masks"]),
    }
    aux["losses"].update(losses)

#####################################################################################

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

def get_ptm(inputs, outputs, interface=False):
  pae = {"residue_weights":inputs["seq_mask"],
         **outputs["predicted_aligned_error"]}
  if interface:
    if "asym_id" not in pae:
      pae["asym_id"] = inputs["asym_id"]
  else:
    if "asym_id" in pae:
      pae.pop("asym_id")
  return confidence.predicted_tm_score(**pae, use_jnp=True)
  
def get_dgram_bins(outputs):
  dgram = outputs["distogram"]["logits"]
  if dgram.shape[-1] == 64:
    dgram_bins = jnp.append(0,jnp.linspace(2.3125,21.6875,63))
  if dgram.shape[-1] == 39:
    dgram_bins = jnp.linspace(3.25,50.75,39) + 1.25
  return dgram_bins

def get_contact_map(outputs, dist=8.0):
  '''get contact map from distogram'''
  dist_logits = outputs["distogram"]["logits"]
  dist_bins = get_dgram_bins(outputs)
  return (jax.nn.softmax(dist_logits) * (dist_bins < dist)).sum(-1)

####################
# confidence metrics
####################
def mask_loss(x, mask=None, mask_grad=False):
  if mask is None:
    return x.mean()
  else:
    x_masked = (x * mask).sum() / (1e-8 + mask.sum())
    if mask_grad:
      return jax.lax.stop_gradient(x.mean() - x_masked) + x_masked
    else:
      return x_masked

def get_exp_res_loss(outputs, mask_1d=None):
  p = jax.nn.sigmoid(outputs["experimentally_resolved"]["logits"])
  p = 1 - p[...,residue_constants.atom_order["CA"]]
  return mask_loss(p, mask_1d)

def get_plddt_loss(outputs, mask_1d=None):
  p = 1 - get_plddt(outputs)
  return mask_loss(p, mask_1d)

def get_pae_loss(outputs, mask_1d=None, mask_1b=None, mask_2d=None):
  p = get_pae(outputs) / 31.0
  p = (p + p.T) / 2
  L = p.shape[0]
  if mask_1d is None: mask_1d = jnp.ones(L)
  if mask_1b is None: mask_1b = jnp.ones(L)
  if mask_2d is None: mask_2d = jnp.ones((L,L))
  mask_2d = mask_2d * mask_1d[:,None] * mask_1b[None,:]
  return mask_loss(p, mask_2d)

def get_con_loss(inputs, outputs, con_opt,
                 mask_1d=None, mask_1b=None, mask_2d=None):

  # get top k
  def min_k(x, k=1, mask=None):
    y = jnp.sort(x if mask is None else jnp.where(mask,x,jnp.nan))
    k_mask = jnp.logical_and(jnp.arange(y.shape[-1]) < k, jnp.isnan(y) == False)
    return jnp.where(k_mask,y,0).sum(-1) / (k_mask.sum(-1) + 1e-8)
  
  # decide on what offset to use
  if "offset" in inputs:
    offset = inputs["offset"]
  else:
    idx = inputs["residue_index"].flatten()
    offset = idx[:,None] - idx[None,:]

  # define distogram
  dgram = outputs["distogram"]["logits"]
  dgram_bins = get_dgram_bins(outputs)

  p = _get_con_loss(dgram, dgram_bins, cutoff=con_opt["cutoff"], binary=con_opt["binary"])
  if "seqsep" in con_opt:
    m = jnp.abs(offset) >= con_opt["seqsep"]
  else:
    m = jnp.ones_like(offset)

  # mask results
  if mask_1d is None: mask_1d = jnp.ones(m.shape[0])
  if mask_1b is None: mask_1b = jnp.ones(m.shape[0])
  
  if mask_2d is None:
    m = jnp.logical_and(m, mask_1b)
  else:
    m = jnp.logical_and(m, mask_2d)  

  p = min_k(p, con_opt["num"], m)
  return min_k(p, con_opt["num_pos"], mask_1d)

def _get_con_loss(dgram, dgram_bins, cutoff=None, binary=True):
  '''dgram to contacts'''
  if cutoff is None: cutoff = dgram_bins[-1]
  bins = dgram_bins < cutoff  
  px = jax.nn.softmax(dgram)
  px_ = jax.nn.softmax(dgram - 1e7 * (1-bins))        
  # binary/cateogorical cross-entropy
  con_loss_cat_ent = -(px_ * jax.nn.log_softmax(dgram)).sum(-1)
  con_loss_bin_ent = -jnp.log((bins * px + 1e-8).sum(-1))
  return jnp.where(binary, con_loss_bin_ent, con_loss_cat_ent)

def get_helix_loss(inputs, outputs, mask_2d=None):  
  # decide on what offset to use
  if "offset" in inputs:
    offset = inputs["offset"]
  else:
    idx = inputs["residue_index"].flatten()
    offset = idx[:,None] - idx[None,:]

  # define distogram
  dgram = outputs["distogram"]["logits"]
  dgram_bins = get_dgram_bins(outputs)

  return _get_helix_loss(dgram, dgram_bins, offset, mask_2d=mask_2d)

def _get_helix_loss(dgram, dgram_bins, offset=None, mask_2d=None, **kwargs):
  '''helix bias loss'''
  x = _get_con_loss(dgram, dgram_bins, cutoff=6.0, binary=True)
  if offset is None:
    if mask_2d is None:
      return jnp.diagonal(x,3).mean()
    else:
      return jnp.diagonal(x * mask_2d,3).sum() + (jnp.diagonal(mask_2d,3).sum() + 1e-8)
  else:
    mask = offset == 3
    if mask_2d is not None:
      mask = jnp.where(mask_2d,mask,0)
    return jnp.where(mask,x,0.0).sum() / (mask.sum() + 1e-8)

####################
# loss functions
####################
def get_dgram_loss(inputs, outputs, copies=1, aatype=None, return_mtx=False):

  batch = inputs["batch"]
  # gather features
  if aatype is None: aatype = batch["aatype"]
  pred = outputs["distogram"]["logits"]

  # get true features
  x, weights = model.modules.pseudo_beta_fn(aatype=aatype,
                                            all_atom_positions=batch["all_atom_positions"],
                                            all_atom_mask=batch["all_atom_mask"])

  dm = jnp.square(x[:,None]-x[None,:]).sum(-1,keepdims=True)
  bin_edges = jnp.linspace(2.3125, 21.6875, pred.shape[-1] - 1)
  true = jax.nn.one_hot((dm > jnp.square(bin_edges)).sum(-1), pred.shape[-1])

  def loss_fn(t,p,m):
    cce = -(t*jax.nn.log_softmax(p)).sum(-1)
    return cce, (cce*m).sum((-1,-2))/(m.sum((-1,-2))+1e-8)
  
  weights = jnp.where(inputs["seq_mask"],weights,0)
  weights_2d = weights[:,None] * weights[None,:]
  weights_2d = jnp.where(inputs["supervised_2d"], weights_2d, 0)
  return _get_pw_loss(true, pred, loss_fn, weights_2d=weights_2d, copies=copies, return_mtx=return_mtx)

def get_fape_loss(inputs, outputs, copies=1, clamp=10.0, return_mtx=False):

  def robust_norm(x, axis=-1, keepdims=False, eps=1e-8):
    return jnp.sqrt(jnp.square(x).sum(axis=axis, keepdims=keepdims) + eps)

  def get_R(N, CA, C):
    (v1,v2) = (C-CA, N-CA)
    e1 = v1 / robust_norm(v1, axis=-1, keepdims=True)
    c = jnp.einsum('li, li -> l', e1, v2)[:,None]
    e2 = v2 - c * e1
    e2 = e2 / robust_norm(e2, axis=-1, keepdims=True)
    e3 = jnp.cross(e1, e2, axis=-1)
    return jnp.concatenate([e1[:,:,None], e2[:,:,None], e3[:,:,None]], axis=-1)

  def get_ij(R,T):
    return jnp.einsum('rji,rsj->rsi',R,T[None,:]-T[:,None])

  def loss_fn(t,p,m):
    fape = robust_norm(t-p)
    fape = jnp.clip(fape, 0, clamp) / 10.0
    return fape, (fape*m).sum((-1,-2))/(m.sum((-1,-2)) + 1e-8)

  true = inputs["batch"]["all_atom_positions"]
  pred = outputs["structure_module"]["final_atom_positions"]

  N,CA,C = (residue_constants.atom_order[k] for k in ["N","CA","C"])

  true_mask = jnp.where(inputs["seq_mask"][:,None],inputs["batch"]["all_atom_mask"],0)
  weights = true_mask[:,N] * true_mask[:,CA] * true_mask[:,C]

  true = get_ij(get_R(true[:,N],true[:,CA],true[:,C]),true[:,CA])
  pred = get_ij(get_R(pred[:,N],pred[:,CA],pred[:,C]),pred[:,CA])

  weights_2d = weights[:,None] * weights[None,:]
  weights_2d = jnp.where(inputs["supervised_2d"], weights_2d, 0)
  return _get_pw_loss(true, pred, loss_fn, weights_2d=weights_2d, copies=copies, return_mtx=return_mtx)

def _get_pw_loss(true, pred, loss_fn, weights=None, copies=1, return_mtx=False, weights_2d=None):
  length = true.shape[0]
  
  if weights is None:
    weights = jnp.ones(length)
  if weights_2d is None:
    weights_2d = weights[:,None] * weights[None,:]
  
  F = {"t":true, "p":pred, "m":weights_2d}  

  if copies > 1:
    (L,C) = (length//copies, copies-1)

    # intra (L,L,F)
    intra = jax.tree_map(lambda x:x[:L,:L], F)
    mtx, loss = loss_fn(**intra)

    # inter (C*L,L,F)
    inter = jax.tree_map(lambda x:x[L:,:L], F)
    if C == 0:
      i_mtx, i_loss = loss_fn(**inter)

    else:
      # (C,L,L,F)
      inter = jax.tree_map(lambda x:x.reshape(C,L,L,-1), inter)
      inter = {"t":inter["t"][:,None],        # (C,1,L,L,F)
               "p":inter["p"][None,:],        # (1,C,L,L,F)
               "m":inter["m"][:,None,:,:,0]}  # (C,1,L,L)             
      
      # (C,C,L,L,F) → (C,C,L,L) → (C,C) → (C) → ()
      i_mtx, i_loss = loss_fn(**inter)
      i_loss = sum([i_loss.min(i).sum() for i in [0,1]]) / 2

    total_loss = (loss + i_loss) / copies
    return (mtx, i_mtx) if return_mtx else total_loss

  else:
    mtx, loss = loss_fn(**F)
    return mtx if return_mtx else loss
  
def get_rmsd_loss(inputs, outputs, L=None, include_L=True, copies=1):
  batch = inputs["batch"]
  true = batch["all_atom_positions"][:,1]
  pred = outputs["structure_module"]["final_atom_positions"][:,1]
  weights = jnp.where(inputs["seq_mask"],batch["all_atom_mask"][:,1],0)
  weights = jnp.where(inputs["interchain_mask"][0],weights,0)
  return _get_rmsd_loss(true, pred, weights=weights, L=L, include_L=include_L, copies=copies)

def _get_rmsd_loss(true, pred, weights=None, L=None, include_L=True, copies=1, eps=1e-8):
  '''
  get rmsd + alignment function
  align based on the first L positions, computed weighted rmsd using all 
  positions (if include_L=True) or remaining positions (if include_L=False).
  '''
  # normalize weights
  length = true.shape[-2]
  if weights is None:
    weights = (jnp.ones(length)/length)[...,None]
  else:
    weights = (weights/(weights.sum(-1,keepdims=True) + eps))[...,None]

  # determine alignment [L]ength and remaining [l]ength
  if copies > 1:
    if L is None:
      L = iL = length // copies; C = copies-1
    else:
      (iL,C) = ((length-L) // copies, copies)
  else:
    (L,iL,C) = (length,0,0) if L is None else (L,length-L,1)

  # slice inputs
  if iL == 0:
    (T,P,W) = (true,pred,weights)
  else:
    (T,P,W) = (x[...,:L,:] for x in (true,pred,weights))
    (iT,iP,iW) = (x[...,L:,:] for x in (true,pred,weights))

  # get alignment and rmsd functions
  (T_mu,P_mu) = ((x*W).sum(-2,keepdims=True)/(W.sum((-1,-2)) + eps) for x in (T,P))
  aln_mtx = _np_kabsch((P-P_mu)*W,(T-T_mu)*W)
  aln_fn = lambda x: (x - P_mu) @ jax.lax.stop_gradient(aln_mtx) + T_mu
  msd_fn = lambda t,p,w: (w*jnp.square(aln_fn(p)-t)).sum(-1).sum(-1)
    
  # compute rmsd
  if iL == 0:
    msd = msd_fn(true,pred,weights)
  elif C > 1:
    # all vs all alignment of remaining, get min RMSD
    iT = iT.reshape(-1,C,1,iL,3).swapaxes(0,-3)
    iP = iP.reshape(-1,1,C,iL,3).swapaxes(0,-3)
    imsd = msd_fn(iT, iP, iW.reshape(-1,C,1,iL,1).swapaxes(0,-3))
    imsd = (imsd.min(0).sum(0) + imsd.min(1).sum(0)) / 2 
    imsd = imsd.reshape(jnp.broadcast_shapes(true.shape[:-2],pred.shape[:-2]))
    msd = (imsd + msd_fn(T,P,W)) if include_L else (imsd/(iW.sum((-1,-2) + eps)))
  else:
    msd = msd_fn(true,pred,weights) if include_L else (msd_fn(iT,iP,iW)/(iW.sum((-1,-2)) + eps))
  rmsd = jnp.sqrt(msd + eps)

  return {"rmsd":rmsd, "align":aln_fn}

def get_sc_loss(inputs, outputs, cfg, mask_1d=None):
  
  batch = inputs["batch"]
  aatype = batch["aatype"] 
  mask = batch["all_atom_mask"]
  if mask_1d is not None:
    mask = jnp.where(mask_1d[:,None],mask,0)

  # atom37 representation
  pred_pos = geometry.Vec3Array.from_array(outputs["structure_module"]["final_atom_positions"])
  pos = geometry.Vec3Array.from_array(batch["all_atom_positions"])

  # atom14 representation
  pred_pos14 = all_atom_multimer.atom37_to_atom14(aatype, pred_pos, mask)[0]
  pos14, mask14, alt_naming_is_better = folding_multimer.compute_atom14_gt(aatype, pos, mask, pred_pos14)

  # frame representation
  pred_frames = all_atom_multimer.atom37_to_frames(aatype, pred_pos, mask)["rigidgroups_gt_frames"]
  frames, frames_mask = folding_multimer.compute_frames(aatype, pos, mask, alt_naming_is_better)
  
  # chi representation
  chi_angles, chi_mask = all_atom_multimer.compute_chi_angles(pos, mask, aatype)
  chi_angles = folding_multimer.get_renamed_chi_angles(aatype, chi_angles, alt_naming_is_better)
  
  # sc_chi
  chi_loss, chi_norm_loss = folding_multimer.supervised_chi_loss(
    sequence_mask=inputs["seq_mask"],
    target_chi_mask=chi_mask,
    target_chi_angles=chi_angles,
    aatype=aatype,
    
    pred_angles=outputs["structure_module"]["sidechains"]['angles_sin_cos'],
    unnormed_angles=outputs["structure_module"]["sidechains"]['unnormalized_angles_sin_cos'],
    
    config=cfg.model.heads.structure_module)[1:]

  # sc_fape
  sc_fape = folding_multimer.sidechain_loss(
    gt_frames=frames,
    gt_frames_mask=frames_mask,
    gt_positions=pos14,
    gt_mask=mask14,

    pred_frames=jax.tree_map(lambda x:x[None], pred_frames),
    pred_positions=jax.tree_map(lambda x:x[None], pred_pos14),
    
    config=cfg.model.heads.structure_module)["loss"]
  
  return {"sc_fape":sc_fape,
          "sc_chi":chi_loss,
          "sc_chi_norm": chi_norm_loss}