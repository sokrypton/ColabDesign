import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.shared.utils import Key
from colabdesign.shared.protein import jnp_rmsd_w, _np_kabsch, _np_rmsd, _np_get_6D_loss
from colabdesign.af.alphafold.model import model, folding, all_atom
from colabdesign.af.alphafold.common import confidence_jax

####################################################
# AF_LOSS - setup loss function
####################################################
class _af_loss:
  # protocol specific loss functions
  def _loss_hallucination(self, inputs, outputs, opt, aux):
    plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
    plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)
    aux["losses"]["plddt"] = plddt_loss.mean()
    self._get_pairwise_loss(inputs, outputs, opt, aux)
    
  def _loss_fixbb(self, inputs, outputs, opt, aux):
    '''get losses'''
    plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
    plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)
    self._get_pairwise_loss(inputs, outputs, opt, aux)

    copies = self._args["copies"]
    if self._args["repeat"] or not self._args["homooligomer"]: copies = 1      
    
    # rmsd loss
    aln = get_rmsd_loss(inputs, outputs, copies=copies)
    rmsd, aux["atom_positions"] = aln["rmsd"], aln["align"](aux["atom_positions"])

    # dgram loss
    aatype = inputs["aatype"]
    dgram_cce = get_dgram_loss(inputs, outputs, copies=copies, aatype=aatype)
    
    aux["losses"].update({"fape": get_fape_loss(inputs, outputs, model_config=self._cfg),
                          "rmsd": rmsd, "dgram_cce": dgram_cce, "plddt":plddt_loss.mean()})

  def _loss_binder(self, inputs, outputs, opt, aux):
    '''get losses'''
    plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
    plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)
    aux["losses"]["plddt"] = plddt_loss[...,self._target_len:].mean()
    self._get_pairwise_loss(inputs, outputs, opt, aux, interface=True)

    if self._args["redesign"]:      
      aln = get_rmsd_loss(inputs, outputs, L=self._target_len, include_L=False)
      align_fn = aln["align"]

      fape = get_fape_loss(inputs, outputs, model_config=self._cfg)
      
      # compute cce of binder + interface
      aatype = inputs["aatype"]
      cce = get_dgram_loss(inputs, outputs, aatype=aatype, return_cce=True)

      aux["losses"].update({"rmsd":aln["rmsd"], "fape":fape, "dgram_cce":cce[self._target_len:,:].mean()})
    
    else:
      align_fn = get_rmsd_loss(inputs, outputs, L=self._target_len)["align"]

    aux["atom_positions"] = align_fn(aux["atom_positions"])

  def _loss_partial(self, inputs, outputs, opt, aux):
    '''get losses'''    
    batch = inputs["batch"]
    plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
    plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)
    aux["losses"]["plddt"] = plddt_loss.mean()
    self._get_pairwise_loss(inputs, outputs, opt, aux)

    def sub(x, p, axis=0):
      fn = lambda y:jnp.take(y,p,axis)
      # fn = lambda y:jnp.tensordot(p,y,(-1,axis)).swapaxes(axis,0)
      return jax.tree_map(fn, x)

    pos = opt["pos"]
    aatype = inputs["aatype"]
    _config = self._cfg.model.heads.structure_module

    # dgram
    dgram = sub(sub(outputs["distogram"]["logits"],pos),pos,1)
    if aatype is not None: aatype = sub(aatype,pos,0)
    aux["losses"]["dgram_cce"] = get_dgram_loss(inputs, pred=dgram, copies=1, aatype=aatype)

    # rmsd
    true = batch["all_atom_positions"]
    pred = sub(outputs["structure_module"]["final_atom_positions"],pos)
    aln = _get_rmsd_loss(true[:,1], pred[:,1])
    aux["losses"]["rmsd"] = aln["rmsd"]

    # fape
    fape_loss = {"loss":0.0}
    struct = outputs["structure_module"]
    traj = {"traj":sub(struct["traj"],pos,-2)}
    folding.backbone_loss(fape_loss, batch, traj, _config)
    aux["losses"]["fape"] = fape_loss["loss"]

    # sidechain specific losses
    if self._args["use_sidechains"]:
      # sc_fape
      pred_pos = sub(struct["final_atom14_positions"],pos)
      sc_struct = {**folding.compute_renamed_ground_truth(self._sc["batch"], pred_pos),
                   "sidechains":{k: sub(struct["sidechains"][k],pos,1) for k in ["frames","atom_pos"]}}

      aux["losses"]["sc_fape"] = folding.sidechain_loss(batch, sc_struct, _config)["loss"]

      # sc_rmsd
      true_pos = all_atom.atom37_to_atom14(batch["all_atom_positions"], self._sc["batch"])
      aln = get_sc_rmsd(true_pos, pred_pos, self._sc["pos"])
      aux["losses"]["sc_rmsd"] = aln["rmsd"]

    # align final atoms
    aux["atom_positions"] = aln["align"](aux["atom_positions"])

  def _get_pairwise_loss(self, inputs, outputs, opt, aux, interface=False):
    '''get pairwise loss features'''

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
    if not interface:
      aux["losses"].update({"con":get_con_loss(dgram, dgram_bins, offset=offset, **opt["con"]).mean(),
                            "helix":get_helix_loss(dgram, dgram_bins, offset=offset, **opt["con"]),
                            "pae":pae.mean()})
    else:
      # split pae/con into inter/intra
      if self.protocol == "binder":
        (L,H) = (self._target_len, opt.get("pos",None))
      else:
        (L,H) = (self._len, None)
      
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
        if self.protocol == "binder":
          return bb,abba.swapaxes(0,1)
        else:
          return aa,abba

      x_offset, ix_offset = split_feats(jnp.abs(offset))
      for k,v in zip(["pae","con"], [pae,dgram]):
        x, ix = split_feats(v)        
        if k == "con":
          aux["losses"]["helix"] = get_helix_loss(x, dgram_bins, x_offset)
          x = get_con_loss(x, dgram_bins, offset=x_offset, **opt["con"])
          ix = get_con_loss(ix, dgram_bins, offset=ix_offset, **opt["i_con"])

        aux["losses"].update({k:x.mean(),f"i_{k}":ix.mean()})

#####################################################################################
#####################################################################################

def get_pw_con_loss(dgram, dgram_bins, cutoff=None, binary=True):
  '''dgram to contacts'''
  if cutoff is None: cutoff = dgram_bins[-1]
  bins = dgram_bins < cutoff  
  px = jax.nn.softmax(dgram)
  px_ = jax.nn.softmax(dgram - 1e7 * (1-bins))        
  # binary/cateogorical cross-entropy
  con_loss_cat_ent = -(px_ * jax.nn.log_softmax(dgram)).sum(-1)
  con_loss_bin_ent = -jnp.log((bins * px + 1e-8).sum(-1))
  return jnp.where(binary, con_loss_bin_ent, con_loss_cat_ent)

def get_con_loss(dgram, dgram_bins, cutoff=None, binary=True,
                 num=1, seqsep=0, offset=None):
  '''convert distogram into contact loss'''  
  x = get_pw_con_loss(dgram, dgram_bins, cutoff, binary)  
  a,b = x.shape
  if offset is None:
    mask = jnp.abs(jnp.arange(a)[:,None] - jnp.arange(b)[None,:]) >= seqsep
  else:
    mask = jnp.abs(offset) >= seqsep
  x = jnp.sort(jnp.where(mask,x,jnp.nan))
  k_mask = (jnp.arange(b) < num) * (jnp.isnan(x) == False)    
  return jnp.where(k_mask,x,0.0).sum(-1) / (k_mask.sum(-1) + 1e-8)

def get_helix_loss(dgram, dgram_bins, offset=None, **kwargs):
  '''helix bias loss'''
  x = get_pw_con_loss(dgram, dgram_bins, cutoff=6.0, binary=True)
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

def get_ptm(outputs):
  pae = outputs["predicted_aligned_error"]
  return confidence_jax.predicted_tm_score_jax(**pae)

####################
# loss functions
####################
def get_dgram_loss(inputs, outputs=None, copies=1, aatype=None, pred=None, return_cce=False):

  batch = inputs["batch"]
  # gather features
  if aatype is None: aatype = batch["aatype"]
  if pred is None: pred = outputs["distogram"]["logits"]

  # get true features
  x, weights = model.modules.pseudo_beta_fn(aatype=aatype,
                                            all_atom_positions=batch["all_atom_positions"],
                                            all_atom_masks=batch["all_atom_mask"])

  dm = jnp.square(x[:,None]-x[None,:]).sum(-1,keepdims=True)
  bin_edges = jnp.linspace(2.3125, 21.6875, pred.shape[-1] - 1)
  true = jax.nn.one_hot((dm > jnp.square(bin_edges)).sum(-1), pred.shape[-1])

  return _get_dgram_loss(true, pred, weights, copies, return_cce=return_cce)

def _get_dgram_loss(true, pred, weights=None, copies=1, return_cce=False):
  
  length = true.shape[0]
  if weights is None: weights = jnp.ones(length)
  F = {"t":true, "p":pred, "m":weights[:,None] * weights[None,:]}  

  def cce_fn(t,p,m):
    cce = -(t*jax.nn.log_softmax(p)).sum(-1)
    return cce, (cce*m).sum((-1,-2))/(m.sum((-1,-2))+1e-8)

  if copies > 1:
    (L,C) = (length//copies, copies-1)

    # intra (L,L,F)
    intra = jax.tree_map(lambda x:x[:L,:L], F)
    cce, cce_loss = cce_fn(**intra)

    # inter (C*L,L,F)
    inter = jax.tree_map(lambda x:x[L:,:L], F)
    if C == 0:
      i_cce, i_cce_loss = cce_fn(**inter)

    else:
      # (C,L,L,F)
      inter = jax.tree_map(lambda x:x.reshape(C,L,L,-1), inter)
      inter = {"t":inter["t"][:,None],        # (C,1,L,L,F)
               "p":inter["p"][None,:],        # (1,C,L,L,F)
               "m":inter["m"][:,None,:,:,0]}  # (C,1,L,L)             
      
      # (C,C,L,L,F) → (C,C,L,L) → (C,C) → (C) → ()
      i_cce, i_cce_loss = cce_fn(**inter)
      i_cce_loss = sum([i_cce_loss.min(i).sum() for i in [0,1]]) / 2

    total_loss = (cce_loss + i_cce_loss) / copies
    return (cce, i_cce) if return_cce else total_loss

  else:
    cce, cce_loss = cce_fn(**F)
    return cce if return_cce else cce_loss

def get_rmsd_loss(inputs, outputs, L=None, include_L=True, copies=1):
  batch = inputs["batch"]
  true = batch["all_atom_positions"][:,1,:]
  pred = outputs["structure_module"]["final_atom_positions"][:,1,:]
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

def get_sc_rmsd(true, pred, sc):
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

#--------------------------------------
# TODO (make copies friendly)
#--------------------------------------
def get_fape_loss(inputs, outputs, model_config, use_clamped_fape=False):
  batch = inputs["batch"]
  sub_batch = jax.tree_map(lambda x: x, batch)
  sub_batch["use_clamped_fape"] = use_clamped_fape
  loss = {"loss":0.0}    
  _config = model_config.model.heads.structure_module
  folding.backbone_loss(loss, sub_batch, outputs["structure_module"], _config)
  return loss["loss"]

def get_6D_loss(inputs, outputs, **kwargs):
  batch = inputs["batch"]
  true = batch["all_atom_positions"]
  pred = outputs["structure_module"]["final_atom_positions"]
  mask = batch["all_atom_mask"]
  return _np_get_6D_loss(true, pred, mask, **kwargs)