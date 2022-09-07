import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.shared.utils import copy_dict
from colabdesign.af.alphafold.model import config

class _af_crop:
  def _crop(self):
    ''' determine positions to crop '''
    (L, max_L, mode) = (sum(self._lengths), self._args["crop_len"], self._args["crop_mode"])
    
    if max_L is None or max_L >= L: crop = False
    elif self._args["copies"] > 1 and not self._args["repeat"]: crop = False
    elif self.protocol in ["partial","binder"]: crop = False
    elif mode == "dist" and not hasattr(self,"_dist"): crop = False
    else: crop = True
    
    if crop:
      if self.protocol == "fixbb":
        self._tmp["cmap"] = self._dist < self.opt["cmap_cutoff"]
    
      if mode == "slide":
        i = jax.random.randint(self.key(),[],0,(L-max_L)+1)
        p = np.arange(i,i+max_L)      

      if mode == "roll":
        i = jax.random.randint(self.key(),[],0,L)
        p = np.sort(np.roll(np.arange(L),L-i)[:max_L])

      if mode == "dist":
        i = jax.random.randint(self.key(),[],0,(L-max_L)+1)
        p = np.sort(self._dist[i].argsort()[1:][:max_L])

      if mode == "pair":
        # pick random pair of interactig crops
        max_L = max_L // 2

        # pick first crop
        i_range = np.append(np.arange(0,(L-2*max_L)+1),np.arange(max_L,(L-max_L)+1))
        i = jax.random.choice(self.key(),i_range,[])
        
        # pick second crop
        j_range = np.append(np.arange(0,(i-max_L)+1),np.arange(i+max_L,(L-max_L)+1))
        if "cmap" in self._tmp:
          # if contact map defined, bias to interacting pairs
          w = np.array([self._tmp["cmap"][i:i+max_L,j:j+max_L].sum() for j in j_range]) + 1e-8
          j = jax.random.choice(self.key(), j_range, [], p=w/w.sum())
        else:
          j = jax.random.choice(self.key(), j_range, [])
             
        p = np.sort(np.append(np.arange(i,i+max_L),np.arange(j,j+max_L)))

      def callback(self):
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
    
    else:
      callback = None
      p = np.arange(sum(self._lengths))

    self.opt["crop_pos"] = p
    return callback	

def crop_feat(feat, pos):  
  '''
  crop features to specified [pos]itions
  '''
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
    if k == "batch":
      new_feat[k] = crop_feat(feat[k], pos)
    if k in idx:
      for i in idx[k]: new_feat[k] = jnp.take(new_feat[k], pos, i)
  
  return new_feat