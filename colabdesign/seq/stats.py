import jax
import jax.numpy as jnp
import numpy as np

def get_stats(X, X_weight=None, labels=None, add_f_ij=True, add_mf_ij=False, add_c=False):
  '''compute f_i/f_ij/f_ijk given msa '''
  n = None
  if X_weight is None:
    Xn = Xs = X
  else:
    Xn, Xs = X*X_weight[:,n,n], X*jnp.sqrt(X_weight[:,n,n])
  f_i = Xn.sum(0)
  o = {"f_i": f_i / f_i.sum(1,keepdims=True)}
  
  if add_f_ij:
    f_ij = jnp.tensordot(Xs,Xs,[0,0])
    o["f_ij"] = f_ij / f_ij.sum((1,3),keepdims=True)
    if add_c: o["c_ij"] = o["f_ij"] - o["f_i"][:,:,n,n] * o["f_i"][n,n,:,:]

  if labels is not None:
    # compute mixture stats
    if jnp.issubdtype(labels, jnp.integer):
      labels = jax.nn.one_hot(labels,labels.max()+1)
    mf_i = jnp.einsum("nc,nia->cia", labels, Xn)
    o["mf_i"] = mf_i/mf_i.sum((0,2),keepdims=True)
    if add_mf_ij:
      mf_ij = jnp.einsum("nc,nia,njb->ciajb", labels, Xs, Xs)
      o["mf_ij"] = mf_ij/mf_ij.sum((0,2,4),keepdims=True)   
      if add_c: o["mc_ij"] = o["mf_ij"] - o["mf_i"][:,:,:,n,n] * o["mf_i"][:,n,n,:,:]
  return o
    
def get_r(a,b):
  a = jnp.array(a).flatten()
  b = jnp.array(b).flatten()
  return jnp.corrcoef(a,b)[0,1]

def inv_cov(X, X_weight=None):
  X = jnp.asarray(X)
  N,L,A = X.shape
  if X_weight is None:
    num_points = N
  else:
    X_weight = jnp.asarray(X_weight)
    num_points = X_weight.sum()
  c = get_stats(X, X_weight, add_mf_ij=True, add_c=True)["c_ij"]
  c = c.reshape(L*A,L*A)
  shrink = 4.5/jnp.sqrt(num_points) * jnp.eye(c.shape[0])
  ic = jnp.linalg.inv(c + shrink)
  return ic.reshape(L,A,L,A)

def get_mtx(W):
  W = jnp.asarray(W)
  # l2norm of 20x20 matrices (note: we ignore gaps)
  raw = jnp.sqrt(jnp.sum(np.square(W[:,1:,:,1:]),(1,3)))
  raw = raw.at[jnp.diag_indices_from(raw)].set(0)

  # apc (average product correction)
  ap = raw.sum(0,keepdims=True) * raw.sum(1,keepdims=True) / raw.sum()
  apc = raw - ap
  apc = apc.at[jnp.diag_indices_from(apc)].set(0)
  return raw, apc

def con_auc(true, pred, mask=None):
  '''compute agreement between predicted and measured contact map'''
  true = jnp.asarray(true)
  pred = jnp.asarray(pred)
  if mask is not None:
    mask = jnp.asarray(mask)
    idx = mask.sum(-1) > 0
    true = true[idx,:][:,idx]
    pred = pred[idx,:][:,idx]
  eval_idx = jnp.triu_indices_from(true, 6)
  pred_, true_ = pred[eval_idx], true[eval_idx] 
  L = (jnp.linspace(0.1,1.0,10)*len(true)).astype(jnp.int32)
  sort_idx = jnp.argsort(pred_)[::-1]
  return jnp.asarray([true_[sort_idx[:l]].mean() for l in L])

