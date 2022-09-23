import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.shared.utils import copy_dict
from colabdesign.af.alphafold.model import config

def setup_crops(self, crop_len=32, crop_mode="slide"):

  assert crop_mode in ["slide","roll","pair","dist"]

  # update arguments
  self._args.update({"use_crop":True, "crop_len":crop_len, "crop_mode":crop_mode})

  # definite initial crop_pos
  self._opt["crop_pos"] = self.opt["crop_pos"] = np.arange(crop_len)
  self._opt["full_len"] = self.opt["full_len"] = sum(self._lengths)

  # add distances
  if self.protocol == "fixbb":
    cb_atoms = self._pdb["cb_feat"]["atoms"]
    cb_atoms[self._pdb["cb_feat"]["mask"] == 0,:] = np.nan
    self._dist = np.sqrt(np.square(cb_atoms[:,None] - cb_atoms[None,:]).sum(-1))

  # callback functions
  def pre_callback(inputs, opt):
    '''crop features'''
    p = opt["crop_pos"]
    inputs = crop_feat(inputs, p)

  def post_callback(aux, opt):
    '''uncrop features'''
    p,L = opt["crop_pos"],opt["full_len"]
    p1,p2 = p[:,None],P[None,:]
    
    def uncrop(x, full=0.0, pair=False):
      if pair:
        return jnp.full((L,L) + x.shape[2:], full).at[p1,p2].set(x)
      else:
        return jnp.full((L,) + x.shape[1:], full).at[p].set(x)

    # TODO uncrop all
    aux["cmap"] = uncrop(aux["cmap"], pair=True)
    aux["pae"] =  uncrop(aux["pae"], np.nan, pair=True)
    x = aux["prev"]
    aux["prev"] = {"prev_pos":  uncrop(x["prev_pos"]),
                   "prev_pair": uncrop(x["prev_pair"], pair=True),
                   "prev_msa_first_row": uncrop(x["prev_msa_first_row"])}

  def design_callback(self):
    (L, max_L, mode) = (sum(self._lengths), self._args["crop_len"], self._args["crop_mode"])

    crop = self._args["use_crop"]
    if max_L is None or max_L >= L: crop = False
    elif self._args["copies"] > 1 and not self._args["repeat"]: crop = False
    elif self.protocol in ["partial","binder"]: crop = False
    elif mode == "dist" and not hasattr(self,"_dist"): crop = False

    if crop:
      if self.protocol == "fixbb":
        self._tmp["cmap"] = self._dist < self.opt["cmap_cutoff"]
    
      if mode == "slide":
        i = np.random.randint(0,(L-max_L)+1)
        p = np.arange(i,i+max_L)      

      if mode == "roll":
        i = np.random.randint(0,L)
        p = np.sort(np.roll(np.arange(L),L-i)[:max_L])

      if mode == "dist":
        i = np.random.randint(0,(L-max_L)+1)
        p = np.sort(self._dist[i].argsort()[1:][:max_L])

      if mode == "pair":
        # pick random pair of interactig crops
        max_L = max_L // 2

        # pick first crop
        i_range = np.append(np.arange(0,(L-2*max_L)+1),np.arange(max_L,(L-max_L)+1))
        i = np.random.choice(i_range)
        
        # pick second crop
        j_range = np.append(np.arange(0,(i-max_L)+1),np.arange(i+max_L,(L-max_L)+1))
        if "cmap" in self._tmp:
          # if contact map defined, bias to interacting pairs
          w = np.array([self._tmp["cmap"][i:i+max_L,j:j+max_L].sum() for j in j_range]) + 1e-8
          j = np.random.choice(j_range, p=w/w.sum())
        else:
          j = np.random.choice(j_range)
             
        p = np.sort(np.append(np.arange(i,i+max_L),np.arange(j,j+max_L)))

    else:
      p = np.arange(sum(self._lengths))
    
    self.opt["crop_pos"] = p
    
    # function to apply after run
    cmap, pae = (np.array(self.aux[k]) for k in ["cmap","pae"])
    mask = np.isnan(pae)
    b = 0.9        
    _pae = self._tmp.get("pae",np.full_like(pae, 31.0))
    self._tmp["pae"] = np.where(mask, _pae, (1-b)*pae + b*_pae)
    
    if self.protocol == "hallucination":
      _cmap = self._tmp.get("cmap",np.zeros_like(cmap))
      self._tmp["cmap"] = np.where(mask, _cmap, (1-b)*cmap + b*_cmap)
    
    self.aux.update(self._tmp)

  # populate callbacks
  for k,v in self._callbacks.items():
    if v is None: v = []
    if not isinstance(v,list): v = [v]
    if k == "pre":    v.append(pre_callback)
    if k == "post":   v.append(post_callback)
    if k == "design": v.append(design_callback)
    self._callbacks[k] = v
    
def crop_feat(feat, pos):  
  '''crop features to specified [pos]itions'''
  if feat is None: return None
  def find(x,k):
    i = []
    for j,y in enumerate(x):
      if y == k: i.append(j)
    return i
  shapes = config.CONFIG.data.eval.feat
  NUM_RES = "num residues placeholder"
  idx = {k:find(v,NUM_RES) for k,v in shapes.items()}
  new_feat = copy_dict(feat)
  for k in new_feat.keys():
    if k in ["batch","prev"]:
      new_feat[k] = crop_feat(feat[k], pos)
    if k in idx:
      for i in idx[k]: new_feat[k] = jnp.take(new_feat[k], pos, i)
  return new_feat