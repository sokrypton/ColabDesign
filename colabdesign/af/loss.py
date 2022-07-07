import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.shared.utils import Key
from colabdesign.shared.protein import jnp_rmsd_w, _np_kabsch, _np_rmsd, _np_get_6D_loss
from colabdesign.af.alphafold.model import model, folding, all_atom
from colabdesign.af.alphafold.common import residue_constants

idx_to_resname = dict((v,k) for k,v in residue_constants.resname_to_idx.items())

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

    copies = 1 if self._args["repeat"] else self._copies
    aatype = inputs["aatype"][0]
    aux["losses"].update({"fape": get_fape_loss(self._batch, outputs, model_config=self._config),
                          "rmsd": get_rmsd_loss_w(self._batch, outputs, copies=copies),
                          "dgram_cce": get_dgram_loss(self._batch, outputs, model_config=self._config, aatype=aatype, copies=copies),
                          "plddt":plddt_loss.mean()})

  def _loss_binder(self, inputs, outputs, opt, aux):
    '''get losses'''
    plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
    plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)
    aux["losses"]["plddt"] = plddt_loss[...,self._target_len:].mean()
    self._get_pairwise_loss(inputs, outputs, opt, aux, interface=True)

    if self._args["redesign"]:
      aatype = inputs["aatype"][0]
      true = self._batch["all_atom_positions"][:,1,:]
      pred = outputs["structure_module"]["final_atom_positions"][:,1,:]
      L = self._target_len
      p = true - true[:L].mean(0)
      q = pred - pred[:L].mean(0)
      p = p @ _np_kabsch(p[:L], q[:L])
      rmsd_binder = jnp.sqrt(jnp.square(p[L:]-q[L:]).sum(-1).mean(-1) + 1e-8)
      fape = get_fape_loss(self._batch, outputs, model_config=self._config)
      aux["losses"].update({"rmsd":rmsd_binder, "fape":fape})
      
      # compute cce of binder + interface
      errors = get_dgram_loss(self._batch, outputs, model_config=self._config, aatype=aatype, copies=1, return_errors=True)
      aux["losses"]["dgram_cce"] = errors["errors"][L:,:].mean()       

  def _loss_partial(self, inputs, outputs, opt, aux):
    '''get losses'''    
    plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
    plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)
    aux["losses"]["plddt"] = plddt_loss.mean()
    self._get_pairwise_loss(inputs, outputs, opt, aux)


    def sub(x, p, axis=0):
      fn = lambda y:jnp.take(y,p,axis)
      # fn = lambda y:jnp.tensordot(p,y,(-1,axis)).swapaxes(axis,0)
      return jax.tree_map(fn, x)

    pos = opt["pos"]
    aatype = inputs["aatype"][0]
    _config = self._config.model.heads.structure_module

    # dgram
    dgram = sub(sub(outputs["distogram"]["logits"],pos),pos,1)
    if aatype is not None: aatype = sub(aatype,pos,0)
    aux["losses"]["dgram_cce"] = get_dgram_loss(self._batch, outputs, model_config=self._config,
                                         logits=dgram, aatype=aatype)

    true = self._batch["all_atom_positions"]
    pred = sub(outputs["structure_module"]["final_atom_positions"],pos)

    # rmsd
    aux["losses"]["rmsd"] = _np_rmsd(true[:,1], pred[:,1])

    # fape
    fape_loss = {"loss":0.0}
    struct = outputs["structure_module"]
    traj = {"traj":sub(struct["traj"],pos,-2)}
    folding.backbone_loss(fape_loss, self._batch, traj, _config)
    aux["losses"]["fape"] = fape_loss["loss"]

    # sidechain specific losses
    if self._args["use_sidechains"]:
      # sc_fape
      pred_pos = sub(struct["final_atom14_positions"],pos)
      sc_struct = {**folding.compute_renamed_ground_truth(self._batch, pred_pos),
                   "sidechains":{k: sub(struct["sidechains"][k],pos,1) for k in ["frames","atom_pos"]}}

      aux["losses"]["sc_fape"] = folding.sidechain_loss(self._batch, sc_struct, _config)["loss"]

      # sc_rmsd
      true_aa = self._batch["aatype"]
      true_pos = all_atom.atom37_to_atom14(self._batch["all_atom_positions"],self._batch)
      aux["losses"]["sc_rmsd"] = get_sc_rmsd(true_pos, pred_pos, true_aa,
                                             atoms_to_exclude=self._args["atoms_to_exclude"])    

  def _get_pairwise_loss(self, inputs, outputs, opt, aux, interface=False):
    '''get pairwise loss features'''

    # decide on what offset to use
    if "offset" in inputs:
      offset = inputs["offset"][0]
    else:
      idx = inputs["residue_index"][0]
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
      L = self._target_len if self.protocol == "binder" else self._len
      H = self._hotspot if hasattr(self,"_hotspot") else None
      
      def split_feats(v):
        '''split pairwise features into intra and inter features'''
        if v is None: return None,None
        aa,bb = v[:L,:L], v[L:,L:]
        ab = v[:L,L:] if H is None else v[H,L:]
        ba = v[L:,:L] if H is None else v[L:,H]
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
          # dgram, dgram_bins, cutoff=None, num=1, seqsep=0, binary=True, offset=None
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