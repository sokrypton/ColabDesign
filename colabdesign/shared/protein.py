import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.af.alphafold.common import residue_constants
from string import ascii_uppercase, ascii_lowercase
alphabet_list = list(ascii_uppercase+ascii_lowercase)

MODRES = {'MSE':'MET','MLY':'LYS','FME':'MET','HYP':'PRO',
          'TPO':'THR','CSO':'CYS','SEP':'SER','M3L':'LYS',
          'HSK':'HIS','SAC':'SER','PCA':'GLU','DAL':'ALA',
          'CME':'CYS','CSD':'CYS','OCS':'CYS','DPR':'PRO',
          'B3K':'LYS','ALY':'LYS','YCM':'CYS','MLZ':'LYS',
          '4BF':'TYR','KCX':'LYS','B3E':'GLU','B3D':'ASP',
          'HZP':'PRO','CSX':'CYS','BAL':'ALA','HIC':'HIS',
          'DBZ':'ALA','DCY':'CYS','DVA':'VAL','NLE':'LEU',
          'SMC':'CYS','AGM':'ARG','B3A':'ALA','DAS':'ASP',
          'DLY':'LYS','DSN':'SER','DTH':'THR','GL3':'GLY',
          'HY3':'PRO','LLP':'LYS','MGN':'GLN','MHS':'HIS',
          'TRQ':'TRP','B3Y':'TYR','PHI':'PHE','PTR':'TYR',
          'TYS':'TYR','IAS':'ASP','GPL':'LYS','KYN':'TRP',
          'CSD':'CYS','SEC':'CYS'}

def pdb_to_string(pdb_file, chains=None, models=None):
  '''read pdb file and return as string'''

  if chains is not None:
    if "," in chains: chains = chains.split(",")
    if not isinstance(chains,list): chains = [chains]
  if models is not None:
    if not isinstance(models,list): models = [models]

  modres = {**MODRES}
  lines = []
  seen = []
  model = 1
  
  if "\n" in pdb_file:
    old_lines = pdb_file.split("\n")
  else:
    with open(pdb_file,"rb") as f:
      old_lines = [line.decode("utf-8","ignore").rstrip() for line in f]  
  for line in old_lines:
    if line[:5] == "MODEL":
      model = int(line[5:])
    if models is None or model in models:
      if line[:6] == "MODRES":
        k = line[12:15]
        v = line[24:27]
        if k not in modres and v in residue_constants.restype_3to1:
          modres[k] = v
      if line[:6] == "HETATM":
        k = line[17:20]
        if k in modres:
          line = "ATOM  "+line[6:17]+modres[k]+line[20:]
      if line[:4] == "ATOM":
        chain = line[21:22]
        if chains is None or chain in chains:
          atom = line[12:12+4].strip()
          resi = line[17:17+3]
          resn = line[22:22+5].strip()
          if resn[-1].isalpha(): # alternative atom
            resn = resn[:-1]
            line = line[:26]+" "+line[27:]
          key = f"{model}_{chain}_{resn}_{resi}_{atom}"
          if key not in seen: # skip alternative placements
            lines.append(line)
            seen.append(key)
      if line[:5] == "MODEL" or line[:3] == "TER" or line[:6] == "ENDMDL":
        lines.append(line)
  return "\n".join(lines)

def renum_pdb_str(pdb_str, Ls=None, renum=True, offset=1):
  if Ls is not None:
    L_init = 0
    new_chain = {}
    for L,c in zip(Ls, alphabet_list):
      new_chain.update({i:c for i in range(L_init,L_init+L)})
      L_init += L  

  n,num,pdb_out = 0,offset,[]
  resnum_ = None
  chain_ = None
  new_chain_ = new_chain[0]
  for line in pdb_str.split("\n"):
    if line[:4] == "ATOM":
      chain = line[21:22]
      resnum = int(line[22:22+5])
      if resnum_ is None: resnum_ = resnum
      if chain_ is None: chain_ = chain
      if resnum != resnum_ or chain != chain_:
        num += (resnum - resnum_)  
        n += 1
        resnum_,chain_ = resnum,chain
      if Ls is not None:
        if new_chain[n] != new_chain_:
          num = offset
          new_chain_ = new_chain[n]
      N = num if renum else resnum
      if Ls is None: pdb_out.append("%s%4i%s" % (line[:22],N,line[26:]))
      else: pdb_out.append("%s%s%4i%s" % (line[:21],new_chain[n],N,line[26:]))        
  return "\n".join(pdb_out)

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
  
def _np_get_6D(all_atom_positions, all_atom_mask=None, use_jax=True, for_trrosetta=False):
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
  n0 = N[...,:,None,:]
  a0,a1 = A[...,:,None,:],A[...,None,:,:]
  b0,b1 = B[...,:,None,:],B[...,None,:,:]
  
  if for_trrosetta:
    out.update({"dist":  _np_len(b0,b1,       use_jax=use_jax),
                "phi":   _np_ang(a0,b0,b1,    use_jax=use_jax, use_acos=True),
                "omega": _np_dih(a0,b0,b1,a1, use_jax=use_jax, use_atan2=True),
                "theta": _np_dih(n0,a0,b0,b1, use_jax=use_jax, use_atan2=True)})  
  else:
    out.update({"dist":  _np_len(b0,b1,       use_jax=use_jax),
                "phi":   _np_ang(a0,b0,b1,    use_jax=use_jax, use_acos=False),
                "omega": _np_dih(a0,b0,b1,a1, use_jax=use_jax, use_atan2=False),
                "theta": _np_dih(n0,a0,b0,b1, use_jax=use_jax, use_atan2=False)})  
  return out

####################
# losses
####################

# RMSD
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

# 6D (see TrRosetta paper)
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

def _np_get_6D_binned(all_atom_positions, all_atom_mask, use_jax=None):
  # TODO: make differentiable, add use_jax option
  ref = _np_get_6D(all_atom_positions,
                   all_atom_mask,
                   use_jax=False, for_trrosetta=True)
  ref = jax.tree_map(jnp.squeeze,ref)

  def mtx2bins(x_ref, start, end, nbins, mask):
    bins = np.linspace(start, end, nbins)
    x_true = np.digitize(x_ref, bins).astype(np.uint8)
    x_true = np.where(mask,0,x_true)
    return np.eye(nbins+1)[x_true][...,:-1]

  mask = (ref["dist"] > 20) | (np.eye(ref["dist"].shape[0]) == 1)
  return {"dist": mtx2bins(ref["dist"],    2.0,  20.0,  37,  mask=mask),
          "omega":mtx2bins(ref["omega"], -np.pi, np.pi, 25,  mask=mask),
          "theta":mtx2bins(ref["theta"], -np.pi, np.pi, 25,  mask=mask),
          "phi":  mtx2bins(ref["phi"],      0.0, np.pi, 13,  mask=mask)}