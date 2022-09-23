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

  # add distances
  if self.protocol == "fixbb":
    cb_atoms = self._pdb["cb_feat"]["atoms"]
    cb_atoms[self._pdb["cb_feat"]["mask"] == 0,:] = np.nan
    self._dist = np.sqrt(np.square(cb_atoms[:,None] - cb_atoms[None,:]).sum(-1))

  # function to apply to inputs
  def pre_callback(inputs, opt):
    '''crop features'''

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

    p = opt["crop_pos"]
    inputs.update(crop_feat(inputs, p))

  # function to apply to outputs
  def post_callback(aux, opt):
    '''uncrop features'''
    # TODO uncrop all

    length = sum(self._lengths)

    def uncrop_feat(x, pos, full=0.0, pair=False):
      if pair:
        p1, p2 = pos[:,None], pos[None,:]
        return jnp.full((length,length)+x.shape[2:],full).at[p1,p2].set(x)
      else:
        return jnp.full((length,)+x.shape[1:],full).at[pos].set(x)

    p = opt["crop_pos"]
    x = aux["prev"]
    full_aux = {"cmap":uncrop_feat(aux["cmap"], p, pair=True),
                "pae": uncrop_feat(aux["pae"],  p, full=jnp.nan, pair=True),
                "prev":{"prev_pos":  uncrop_feat(x["prev_pos"],  p),
                        "prev_pair": uncrop_feat(x["prev_pair"], p, pair=True),
                        "prev_msa_first_row": uncrop_feat(x["prev_msa_first_row"], p)},
                }
    if self.protocol == "fixbb":
      full_aux.update({
        "atom_positions":uncrop_feat(aux["atom_positions"],  p),
        "plddt":uncrop_feat(aux["plddt"], p)})

    aux.update(full_aux)

  # function to apply during design
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
    
    if crop:
      # function to apply after run
      cmap, pae, atoms, plddt = (self.aux[k] for k in ["cmap","pae", "atom_positions", "plddt"])
      pae_mask = np.isnan(pae)
      b = 0.5  

      # update pae      
      _pae = self._tmp.get("pae",np.full_like(pae, 31.0))
      self._tmp["pae"] = np.where(pae_mask, _pae, (1-b)*pae + b*_pae)

      if self.protocol == "fixbb" and self._args["realign"]:
        # update atoms
        #if "atom_positions" not in self._tmp:
        #  mu = self._inputs["batch"]["all_atom_positions"]
        #  mu_mask = self._inputs["batch"]["all_atom_mask"][...,None]
        #  mu = (mu * mu_mask).sum(0) / (1e-8 + mu_mask.sum(0))

        #  mu = self._inputs["batch"]["all_atom_positions"] - mu
        #  mu = (mu * mu_mask).sum(0) / (1e-8 + mu_mask.sum(0))

        #  self._tmp["atom_positions"] = np.broadcast_to(mu, atoms.shape)
        #_atoms = self._tmp["atom_positions"]
        #self._tmp["atom_positions"] = np.where(atoms == 0, _atoms, (1-b)*atoms + b*_atoms)

        _plddt = self._tmp.get("plddt",np.zeros_like(plddt))
        self._tmp["plddt"] = np.where(plddt == 0, _plddt, (1-b)*plddt + b*_plddt)
      
      if self.protocol == "hallucination":
        # update contact map
        _cmap = self._tmp.get("cmap",np.zeros_like(cmap))
        self._tmp["cmap"] = np.where(pae_mask, _cmap, (1-b)*cmap + b*_cmap)
      
      self.aux.update(self._tmp)

  # populate callbacks
  for k,v in self._callbacks.items():
    if v is None: v = []
    if not isinstance(v,list): v = [v]
    if k == "pre":    v.append(pre_callback)
    if k == "post":   v.append(post_callback)
    if k == "design": v.append(design_callback)
    self._callbacks[k] = v