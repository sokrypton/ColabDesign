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

    # define masks
    seq_mask_1d = inputs["seq_mask"]
    seq_mask_2d = and_masks(seq_mask_1d[:,None],seq_mask_1d[None,:])
    mask_1d = and_masks(seq_mask_1d,inputs["supervised"],inputs["interchain_mask"][0])
    mask_2d = and_masks(seq_mask_2d,inputs["supervised_2d"],inputs["interchain_mask"])

    fix_mask_1d = inputs["fix_pos"]
    fix_mask_2d = and_masks(fix_mask_1d[:,None],fix_mask_1d[None,:])
    mask_2d = jnp.where(opt["partial_loss"],jnp.where(fix_mask_2d,False,mask_2d),mask_2d)

    # RMSD
    rmsd_flags = dict(include_L=False, L=self._target_len) if self.protocol == "binder" else {}
    aln = get_rmsd_loss(inputs, outputs, copies=a["copies"], mask_1d=mask_1d, **rmsd_flags)
    if a["realign"]: aux["atom_positions"] = aln["align"](aux["atom_positions"]) * aux["atom_mask"][...,None]    

    # add other losses
    aux["losses"].update({
      "fape":      get_fape_loss(inputs, outputs, copies=a["copies"], clamp=opt["fape_cutoff"], mask_2d=mask_2d),
      "dgram_cce": get_dgram_loss(inputs, outputs, copies=a["copies"], aatype=inputs["aatype"], mask_2d=mask_2d),
      "rmsd":      aln["rmsd"],
    })

    if a["use_sidechains"] and "fix_pos" in inputs:
      aux["losses"].update(get_sc_loss(inputs, outputs, mask_1d=inputs["fix_pos"], cfg=self._cfg))    

  def _loss_unsupervised(self, inputs, outputs, aux, mask_1d=None, mask_2d=None):

    opt = inputs["opt"]

    # define masks
    seq_mask_1d = inputs["seq_mask"]
    seq_mask_2d = and_masks(seq_mask_1d[:,None],seq_mask_1d[None,:])
    intra_mask_2d = inputs["asym_id"][:,None] == inputs["asym_id"][None,:]

    # partial masks
    mask_1d = and_masks(seq_mask_1d,inputs["unsupervised"])
    unsup_2d = or_masks(inputs["unsupervised_2d"],inputs["interchain_mask"] == 0)
    mask_2d = and_masks(seq_mask_2d,unsup_2d)
    mask_2d_intra = and_masks(intra_mask_2d,mask_2d)
    mask_2d_inter = and_masks((intra_mask_2d == False),mask_2d)
    masks_partial = {"mask_1d":mask_1d,              "mask_2d":mask_2d,
     "intra_masks": {"mask_1d":mask_2d_intra.any(1), "mask_2d":mask_2d_intra},
     "inter_masks": {"mask_1d":mask_2d_inter.any(1), "mask_2d":mask_2d_inter}}
    
    # full masks
    masks_full =    {"mask_1d":seq_mask_1d, "mask_2d":seq_mask_2d,
      "intra_masks":{"mask_1d":seq_mask_1d, "mask_2d":and_masks(intra_mask_2d,seq_mask_2d)},
      "inter_masks":{"mask_1d":seq_mask_1d, "mask_2d":and_masks((intra_mask_2d == False),seq_mask_2d)}}
    
    masks = jax.tree_map(lambda x,y: jnp.where(opt["partial_loss"],x,y), masks_partial, masks_full)
    if "hotspot" in inputs:
      masks["inter_masks"]["mask_1d"] = and_masks(inputs["hotspot"],seq_mask_1d)

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
def mask_loss(x, mask):
  mask = mask.astype(x.dtype)
  x_masked = (x * mask).sum() / (1e-8 + mask.sum())
  return x_masked

def get_exp_res_loss(outputs, mask_1d):
  p = jax.nn.sigmoid(outputs["experimentally_resolved"]["logits"])
  p = 1 - p[...,residue_constants.atom_order["CA"]]
  return mask_loss(p, mask_1d)

def get_plddt_loss(outputs, mask_1d):
  p = 1 - get_plddt(outputs)
  return mask_loss(p, mask_1d)

def get_pae_loss(outputs, mask_1d=None, mask_1b=None, mask_2d=None):
  p = get_pae(outputs) / 31.0
  p = (p + p.T) / 2
  L = p.shape[0]
  if mask_1d is None: mask_1d = jnp.full(L,True)
  if mask_1b is None: mask_1b = jnp.full(L,True)
  if mask_2d is None: mask_2d = jnp.full((L,L),True)
  mask_2d = and_masks(mask_1d[:,None],mask_1b[None,:],mask_2d)
  return mask_loss(p, mask_2d)

def get_con_loss(inputs, outputs, con_opt,
  mask_1d=None, mask_1b=None, mask_2d=None):

  # get top k
  def min_k(x, k=1, mask=None):
    y = jnp.sort(x if mask is None else jnp.where(mask,x,jnp.nan))
    k_mask = and_masks(jnp.arange(y.shape[-1]) < k, jnp.isnan(y) == False)
    return (k_mask * y).sum(-1) / (k_mask.sum(-1) + 1e-8)
  
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
    m = jnp.full_like(offset, True)

  # mask results
  if mask_2d is not None:
    m = and_masks(m,mask_2d)
  elif mask_1b is not None:
    m = and_masks(m,mask_1b)
  
  if mask_1d is None:
    mask_1d = jnp.full(p.shape[0],True)

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
      return jnp.diagonal(x * mask_2d,3).sum() / (jnp.diagonal(mask_2d,3).sum() + 1e-8)
  else:
    mask = offset == 3
    if mask_2d is not None:
      mask = and_masks(mask, mask_2d)
    return (mask * x).sum() / (mask.sum() + 1e-8)

####################
# loss functions
####################
def get_dgram_loss(inputs, outputs, copies=1, aatype=None, mask_1d=None, mask_2d=None):

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
  
  if mask_1d is not None: weights *= mask_1d
  weights_2d = weights[:,None] * weights[None,:]
  if mask_2d is not None: weights_2d *= mask_2d
  return _get_pw_loss(true, pred, loss_fn, weights_2d, copies=copies)

def get_fape_loss(inputs, outputs, copies=1, clamp=10.0, mask_1d=None, mask_2d=None):

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

  true = get_ij(get_R(true[:,N],true[:,CA],true[:,C]),true[:,CA])
  pred = get_ij(get_R(pred[:,N],pred[:,CA],pred[:,C]),pred[:,CA])

  true_mask = inputs["batch"]["all_atom_mask"]
  weights = true_mask[:,N] * true_mask[:,CA] * true_mask[:,C]
  if mask_1d is not None: weights *= mask_1d
  weights_2d = weights[:,None] * weights[None,:]
  if mask_2d is not None: weights_2d *= mask_2d

  return _get_pw_loss(true, pred, loss_fn, weights_2d, copies=copies)

def _get_pw_loss(true, pred, loss_fn, weights, copies=1):
  dtype = pred.dtype
  length = pred.shape[0]
  F = {"t":true.astype(dtype), "p":pred, "m":weights.astype(dtype)}  

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
    return total_loss

  else:
    mtx, loss = loss_fn(**F)
    return loss
  
def get_rmsd_loss(inputs, outputs, L=None, include_L=True, copies=1, mask_1d=None):
  batch = inputs["batch"]
  true = batch["all_atom_positions"][:,1]
  pred = outputs["structure_module"]["final_atom_positions"][:,1]
  if mask_1d is None: mask_1d = jnp.ones(true.shape[0])
  mask_1d = mask_1d * batch["all_atom_mask"][:,1]
  return _get_rmsd_loss(true, pred, weights=mask_1d, L=L, include_L=include_L, copies=copies)

def _get_rmsd_loss(true, pred, weights, L=None, include_L=True, copies=1, eps=1e-8):
  '''
  get rmsd + alignment function
  align based on the first L positions, computed weighted rmsd using all 
  positions (if include_L=True) or remaining positions (if include_L=False).
  '''
  dtype = pred.dtype
  length = pred.shape[0]  
  true, weights = true.astype(dtype), weights.astype(dtype)
  
  # normalize weights
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

def get_sc_loss(inputs, outputs, mask_1d, cfg):

  batch = inputs["batch"]
  aatype = batch["aatype"] 
  mask = batch["all_atom_mask"] * mask_1d[:,None]

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

def and_masks(*m):
  mask = m[0]
  for n in range(1,len(m)):
    mask = jnp.logical_and(mask,m[n])
  return mask
def or_masks(*m):
  mask = m[0]
  for n in range(1,len(m)):
    mask = jnp.logical_or(mask,m[n])
  return mask

def numpy_callback(x):
  # Need to forward-declare the shape & dtype of the expected output.
  result_shape = jax.core.ShapedArray(x.shape, x.dtype)
  return jax.pure_callback(np.sin, result_shape, x)

def get_chain_indices(chain_boundaries):
    """Returns a list of tuples indicating the start and end indices for each chain."""
    chain_starts_ends = []
    unique_chains = np.unique(chain_boundaries)
    for chain in unique_chains:
        positions = np.where(chain_boundaries == chain)[0]
        chain_starts_ends.append((positions[0], positions[-1]))
    return chain_starts_ends

def get_ifptm(af):
    import string
    from copy import deepcopy
    cmap = get_contact_map(af.aux['debug']['outputs'], 8)  # Define interface with 8A between Cb-s

    # Initialize seq_mask to all False
    inputs_ifptm = {}
    input_pairwise_iptm = {}
    pairwise_iptm = {}

    # Prepare a dictionary to collect results
    pairwise_if_ptm = {}
    chain_starts_ends = get_chain_indices(af._inputs['asym_id'])

    # Generate chain labels (A, B, C, ...)
    chain_labels = list(string.ascii_uppercase)

    for i, (start_i, end_i) in enumerate(chain_starts_ends):
        chain_label_i = chain_labels[i % len(chain_labels)]  # Wrap around if more than 26 chains
        for j, (start_j, end_j) in enumerate(chain_starts_ends):
            chain_label_j = chain_labels[j % len(chain_labels)]  # Wrap around if more than 26 chains
            if i < j:  # Avoid self-comparison and duplicate comparisons
                outputs = deepcopy(af.aux['debug']['outputs'])
                contacts = np.where(cmap[start_i:end_i+1, start_j:end_j+1] >= 0.6)
                key = f"{chain_label_i}-{chain_label_j}"
                
                if contacts[0].size > 0:  # If there are contacts
                    # Convert local chain positions back to global positions using JAX
                    global_i_positions = contacts[0] + start_i
                    global_j_positions = contacts[1] + start_j
                    global_positions = list(set(np.concatenate((global_i_positions, global_j_positions))))
                    global_positions = np.array(global_positions, dtype=int)
                    global_positions.sort()
                    
                    # Initialize new input dictionary
                    inputs_ifptm['seq_mask'] = np.full(len(af._inputs['asym_id']), 0, dtype=float)
                    inputs_ifptm['asym_id'] = af._inputs['asym_id']
                    # Update seq_mask for these positions to True within inputs
                    inputs_ifptm['seq_mask'][global_positions] = 1                
                    # Call get_ptm with updated inputs and outputs
                    pairwise_if_ptm[key] = get_ptm(inputs_ifptm, outputs, interface=True)
                else:
                    pairwise_if_ptm[key] = 0
                
                # Also adding regular i_ptm, pairwise
                outputs = deepcopy(af.aux['debug']['outputs'])
                input_pairwise_iptm['seq_mask'] = np.full(len(af._inputs['asym_id']), 0, dtype=float)
                input_pairwise_iptm['asym_id'] = af._inputs['asym_id']
                input_pairwise_iptm['seq_mask'][np.concatenate((np.arange(start_i,end_i+1), 
                                                               np.arange(start_j,end_j+1)))] = 1  
                pairwise_iptm[key] = get_ptm(input_pairwise_iptm, outputs, interface=True)

    return pairwise_if_ptm, pairwise_iptm
