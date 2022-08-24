import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.shared.utils import Key, copy_dict
from colabdesign.shared.protein import jnp_rmsd_w, _np_kabsch, _np_rmsd, _np_get_6D_loss
from colabdesign.af.alphafold.model import model, folding, all_atom
from colabdesign.af.alphafold.common import confidence_jax, residue_constants

####################################################
# AF_LOSS - setup loss function
####################################################

class _af_loss:
  # protocol specific loss functions
  def _loss_hallucination(self, inputs, outputs, opt, aux):
    aux["losses"]["plddt"] = get_plddt_loss(outputs).mean()
    aux["losses"].update(get_pw_conf_loss(inputs, outputs, opt))
    
  def _loss_fixbb(self, inputs, outputs, opt, aux):
    '''get losses'''
    aux["losses"].update(get_pw_conf_loss(inputs, outputs, opt))

    copies = self._args["copies"] if self._args["homooligomer"] else 1
    
    # rmsd loss
    aln = get_rmsd_loss(inputs, outputs, copies=copies)
    rmsd, aux["atom_positions"] = aln["rmsd"], aln["align"](aux["atom_positions"])

    # dgram loss
    dgram_cce = get_dgram_loss(inputs, outputs, copies=copies, aatype=inputs["aatype"])
    
    # fape loss
    fape = get_fape_loss(inputs, outputs, copies=copies, clamp=opt["fape_cutoff"])
    aux["losses"].update({"rmsd": rmsd, "dgram_cce": dgram_cce, "fape":fape, 
                          "plddt":get_plddt_loss(outputs).mean()})

  def _loss_binder(self, inputs, outputs, opt, aux):
    '''get losses'''
    aux["losses"]["plddt"] = get_plddt_loss(outputs)[...,self._target_len:].mean()
    aux["losses"].update(get_pw_conf_loss(inputs, outputs, opt, L=self._target_len), H=opt.get("pos",None))

    if self._args["redesign"]:      
      aln = get_rmsd_loss(inputs, outputs, L=self._target_len, include_L=False)
      align_fn = aln["align"]
      
      # compute cce of binder + interface
      aatype = inputs["aatype"]
      cce = get_dgram_loss(inputs, outputs, aatype=aatype, return_mtx=True)
      cce = cce[self._target_len:,:].mean()

      # compute fape
      fape = get_fape_loss(inputs, outputs, clamp=opt["fape_cutoff"], return_mtx=True)
      fape = fape[self._target_len:,:].mean()

      aux["losses"].update({"rmsd":aln["rmsd"], "dgram_cce":cce, "fape":fape})
    
    else:
      align_fn = get_rmsd_loss(inputs, outputs, L=self._target_len)["align"]

    aux["atom_positions"] = align_fn(aux["atom_positions"])

  def _loss_partial(self, inputs, outputs, opt, aux):
    '''get losses'''    
    batch = inputs["batch"]
    aux["losses"]["plddt"] = get_plddt_loss(outputs).mean()
    aux["losses"].update(get_pw_conf_loss(inputs, outputs, opt))

    copies = self._args["copies"] if self._args["homooligomer"] else 1

    # subset inputs/outputs
    def sub(x, p, axis=0):
      fn = lambda y:jnp.take(y,p,axis)
      # fn = lambda y:jnp.tensordot(p,y,(-1,axis)).swapaxes(axis,0)
      return jax.tree_map(fn, x)
    
    pos = opt["pos"]
    inputs["aatype"] = sub(inputs["aatype"],pos)
    outputs["structure_module"]["final_atom_positions"] = sub(outputs["structure_module"]["final_atom_positions"],pos)
    outputs["distogram"]["logits"] = sub(sub(outputs["distogram"]["logits"],pos),pos,1)

    # dgram_cce
    aux["losses"]["dgram_cce"] = get_dgram_loss(inputs, outputs, copies=copies, aatype=inputs["aatype"])

    # rmsd
    aln = get_rmsd_loss(inputs, outputs, copies=copies)
    aux["losses"]["rmsd"] = aln["rmsd"]

    # fape
    aux["losses"]["fape"] = get_fape_loss(inputs, outputs, copies=copies, clamp=opt["fape_cutoff"])

    # sidechain specific losses
    if self._args["use_sidechains"] and copies == 1:
    
      _config = self._cfg.model.heads.structure_module
      # sc_fape
      struct = outputs["structure_module"]
      pred_pos = sub(struct["final_atom14_positions"],pos)
      
      if not self._args["use_multimer"]:
        sc_struct = {**folding.compute_renamed_ground_truth(self._sc["batch"], pred_pos),
                     "sidechains":{k: sub(struct["sidechains"][k],pos,1) for k in ["frames","atom_pos"]}}
        aux["losses"]["sc_fape"] = folding.sidechain_loss(batch, sc_struct, _config)["loss"]

      else:
        
        # TODO
        print("ERROR: 'sc_fape' not currently supported for 'multimer' mode")
        aux["losses"]["sc_fape"] = 0.0

      # sc_rmsd
      true_pos = all_atom.atom37_to_atom14(batch["all_atom_positions"], self._sc["batch"])
      aln = _get_sc_rmsd(true_pos, pred_pos, self._sc["pos"])
      aux["losses"]["sc_rmsd"] = aln["rmsd"]

    # align final atoms
    aux["atom_positions"] = aln["align"](aux["atom_positions"])

#####################################################################################

def _get_pw_con_loss(dgram, dgram_bins, cutoff=None, binary=True):
  '''dgram to contacts'''
  if cutoff is None: cutoff = dgram_bins[-1]
  bins = dgram_bins < cutoff  
  px = jax.nn.softmax(dgram)
  px_ = jax.nn.softmax(dgram - 1e7 * (1-bins))        
  # binary/cateogorical cross-entropy
  con_loss_cat_ent = -(px_ * jax.nn.log_softmax(dgram)).sum(-1)
  con_loss_bin_ent = -jnp.log((bins * px + 1e-8).sum(-1))
  return jnp.where(binary, con_loss_bin_ent, con_loss_cat_ent)

def _get_con_loss(dgram, dgram_bins, cutoff=None, binary=True,
                 num=1, seqsep=0, offset=None):
  '''convert distogram into contact loss'''  
  x = _get_pw_con_loss(dgram, dgram_bins, cutoff, binary)  
  a,b = x.shape
  if offset is None:
    mask = jnp.abs(jnp.arange(a)[:,None] - jnp.arange(b)[None,:]) >= seqsep
  else:
    mask = jnp.abs(offset) >= seqsep
  x = jnp.sort(jnp.where(mask,x,jnp.nan))
  k_mask = (jnp.arange(b) < num) * (jnp.isnan(x) == False)    
  return jnp.where(k_mask,x,0.0).sum(-1) / (k_mask.sum(-1) + 1e-8)

def _get_helix_loss(dgram, dgram_bins, offset=None, **kwargs):
  '''helix bias loss'''
  x = _get_pw_con_loss(dgram, dgram_bins, cutoff=6.0, binary=True)
  if offset is None:
    return jnp.diagonal(x,3).mean()
  else:
    mask = offset == 3
    return jnp.where(mask,x,0.0).sum() / (mask.sum() + 1e-8)

def get_contact_map(outputs, dist=8.0):
  '''get contact map from distogram'''
  dist_logits = outputs["distogram"]["logits"]
  dist_bins = jax.numpy.append(0,outputs["distogram"]["bin_edges"])
  dist_mtx = dist_bins[dist_logits.argmax(-1)]
  return (jax.nn.softmax(dist_logits) * (dist_bins < dist)).sum(-1)

####################
# confidence metrics
####################
def get_plddt_loss(outputs):
  p = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
  return (p * jnp.arange(p.shape[-1])[::-1]).mean(-1)

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
  pae = outputs["predicted_aligned_error"]
  if "asym_id" not in pae:
    pae["asym_id"] = inputs["asym_id"]
  return confidence_jax.predicted_tm_score_jax(**pae, interface=interface)

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
  
  return _get_pw_loss(true, pred, loss_fn, weights=weights, copies=copies, return_mtx=return_mtx)

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

  true_mask = inputs["batch"]["all_atom_mask"]
  weights = true_mask[:,N] * true_mask[:,CA] * true_mask[:,C]

  true = get_ij(get_R(true[:,N],true[:,CA],true[:,C]),true[:,CA])
  pred = get_ij(get_R(pred[:,N],pred[:,CA],pred[:,C]),pred[:,CA])

  return _get_pw_loss(true, pred, loss_fn, weights=weights, copies=copies, return_mtx=return_mtx)

def _get_pw_loss(true, pred, loss_fn, weights=None, copies=1, return_mtx=False):
  length = true.shape[0]
  
  if weights is None:
    weights = jnp.ones(length)
  
  F = {"t":true, "p":pred, "m":weights[:,None] * weights[None,:]}  

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
  weights = batch["all_atom_mask"][:,1]
  return _get_rmsd_loss(true, pred, weights=weights, L=L, include_L=include_L, copies=copies)

def _get_rmsd_loss(true, pred, weights=None, L=None, include_L=True, copies=1):
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
    weights = (weights/(weights.sum(-1,keepdims=True) + 1e-8))[...,None]

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
  (T_mu,P_mu) = ((x*W).sum(-2,keepdims=True)/W.sum((-1,-2)) for x in (T,P))
  aln = _np_kabsch((P-P_mu)*W, T-T_mu)   
  align_fn = lambda x: (x - P_mu) @ aln + T_mu
  msd_fn = lambda t,p,w: (w*jnp.square(align_fn(p)-t)).sum((-1,-2))
  
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
    msd = (imsd + msd_fn(T,P,W)) if include_L else (imsd/iW.sum((-1,-2)))
  else:
    msd = msd_fn(true,pred,weights) if include_L else (msd_fn(iT,iP,iW)/iW.sum((-1,-2)))
  rmsd = jnp.sqrt(msd + 1e-8)

  return {"rmsd":rmsd, "align":align_fn}

def _get_sc_rmsd(true, pred, sc):
  '''get sidechain rmsd + alignment function'''

  # select atoms
  (T, P) = (true.reshape(-1,3), pred.reshape(-1,3))
  (T, T_alt, P) = (T[sc["pos"]], T[sc["pos_alt"]], P[sc["pos"]])

  # select non-ambigious atoms
  (T_na, P_na) = (T[sc["non_amb"]], P[sc["non_amb"]])

  # get alignment of non-ambigious atoms
  if "weight_non_amb" in sc:
    T_mu_na = (T_na * sc["weight_non_amb"]).sum(0)
    P_mu_na = (P_na * sc["weight_non_amb"]).sum(0)
    aln = _np_kabsch((P_na-P_mu_na) * sc["weight_non_amb"], T_na-T_mu_na)
  else:
    T_mu_na, P_mu_na = T_na.mean(0), P_na.mean(0)
    aln = _np_kabsch(P_na-P_mu_na, T_na-T_mu_na)

  # apply alignment to all atoms
  align_fn = lambda x: (x - P_mu_na) @ aln + T_mu_na
  P = align_fn(P)

  # compute rmsd
  sd = jnp.minimum(jnp.square(P-T).sum(-1), jnp.square(P-T_alt).sum(-1))
  if "weight" in sc:
    msd = (sd*sc["weight"]).sum()
  else:
    msd = sd.mean()
  rmsd = jnp.sqrt(msd + 1e-8)
  return {"rmsd":rmsd, "align":align_fn}


def get_pw_conf_loss(inputs, outputs, opt, L=None, H=None):
  '''get pairwise confidence loss (pae, contact, helix)'''

  # decide on what offset to use
  if "offset" in inputs:
    offset = inputs["offset"]
  else:
    idx = inputs["residue_index"]
    offset = idx[:,None] - idx[None,:]

  # pae loss
  pae_prob = jax.nn.softmax(outputs["predicted_aligned_error"]["logits"])
  pae = (pae_prob * jnp.arange(pae_prob.shape[-1])).mean(-1)
  
  # define distogram
  dgram = outputs["distogram"]["logits"]
  dgram_bins = jnp.append(0,outputs["distogram"]["bin_edges"])
  
  if L is None:
    return {"con":_get_con_loss(dgram, dgram_bins, offset=offset, **opt["con"]).mean(),
            "helix":_get_helix_loss(dgram, dgram_bins, offset=offset, **opt["con"]),
            "pae":pae.mean()}
  else:
    # split pae/con into inter/intra    
    def split_feats(v):
      '''split pairwise features into intra and inter features'''
      if v is None: 
        return None,None
      
      (aa,bb) =   (v[:L,:L],v[L:,L:])
      if H is None:
        (ab,ba) = (v[:L,L:],v[L:,:L])
      else:
        (ab,ba) = (v[H, L:],v[L:, H])

      abba = (ab + ba.swapaxes(0,1)) / 2
      
      # TODO (aa vs bb)
      return bb, abba.swapaxes(0,1)

    losses = {}
    x_offset, ix_offset = split_feats(jnp.abs(offset))
    for k,v in zip(["pae","con"], [pae,dgram]):
      x, ix = split_feats(v)        
      if k == "con":
        losses["helix"] = _get_helix_loss(x, dgram_bins, x_offset)
        x = _get_con_loss(x, dgram_bins, offset=x_offset, **opt["con"])
        ix = _get_con_loss(ix, dgram_bins, offset=ix_offset, **opt["i_con"])
        
      losses.update({k:x.mean(),
                     f"i_{k}":ix.mean()})
    return losses

def get_seq_ent_loss(inputs, outputs, opt):
  x = inputs["seq"]["logits"] / opt["temp"]
  ent = -(jax.nn.softmax(x) * jax.nn.log_softmax(x)).sum(-1)
  if "pos" in opt and "fix_seq" in opt:
    ent = ent.at[...,opt["pos"]].set(0.0)
  return {"seq_ent":ent.mean()}