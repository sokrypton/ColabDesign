import jax
import jax.numpy as jnp
import numpy
from math import log

def _kmeans(X, X_weight, n_clusters=8, n_init=10, max_iter=300, tol=1e-4, seed=0):
  '''kmeans implemented in jax'''
  
  def _dist(a,b):
    sm = a @ b.T

    a_norm = jnp.square(a).sum(-1)
    b_norm = jnp.square(b).sum(-1)
    
    return jnp.abs(a_norm[:,None] + b_norm[None,:] - 2 * sm)

  def _kmeans_plus_plus(key, X, X_weight, n_clusters):
    '''kmeans++ implemented in jax, for initialization'''
    n_samples, n_features = X.shape
    n_candidates = 2 + int(log(n_clusters))
    
    def loop(m,c):
      n,k = c
      
      inf_mask = jnp.inf * (jnp.arange(n_clusters) > n)
      p = (inf_mask + _dist(X,m)).min(-1)

      # sample candidates
      candidates = jax.random.choice(k, jnp.arange(n_samples),
                                     shape=(n_candidates,),
                                     p=p/p.sum(), replace=False)
      
      # pick sample that decreases inertia the most
      dist = jnp.minimum(p[:,None],_dist(X,X[candidates]))
      i = candidates[(X_weight[:,None] * dist).sum(0).argmin()]
      return m.at[n].set(X[i]), None

    i = jax.random.choice(key,jnp.arange(n_samples))
    init_means = jnp.zeros((n_clusters,n_features)).at[0].set(X[i])
    carry = (jnp.arange(1,n_clusters), jax.random.split(key, n_clusters-1))
    return jax.lax.scan(loop, init_means, carry)[0]

  def _E(means):
    # get labels
    return _dist(X,means).argmin(-1)

  def _M(labels):
    # get means
    labels = jax.nn.one_hot(labels, n_clusters)
    labels = labels * X_weight[:,None]
    labels /= labels.sum(0) + 1e-8
    return labels.T @ X
  
  def _inertia(means):
    # compute score: sum(min(dist(X,means)))
    sco = _dist(X,means).min(-1)
    return (X_weight * sco).sum()

  def single_run(key):
    # initialize
    init_means = _kmeans_plus_plus(key, X, X_weight, n_clusters)

    # run EM
    if tol == 0:
      means = jax.lax.scan(lambda mu,_:(_M(_E(mu)),None), init_means,
                           None, length=max_iter)[0]
    else:
      def EM(x):
        old_mu, old_sco, _, n = x
        new_mu = _M(_E(old_mu))
        new_sco = _inertia(new_mu)
        return new_mu, new_sco, old_sco, n+1
      def check(x):
        _, new_sco, old_sco, n = x
        return ((old_sco-new_sco) > tol) & (n < max_iter)
      init = EM((init_means,jnp.inf,None,0))
      means = jax.lax.while_loop(check, EM, init)[0]

    return {"labels":_E(means),
            "means":means,
            "inertia":_inertia(means)}

  # mulitple runs
  key = jax.random.PRNGKey(seed)
  if n_init > 0:
    out = jax.vmap(single_run)(jax.random.split(key,n_init))
    i = out["inertia"].argmin()
    out = jax.tree_map(lambda x:x[i],out)
  else:
    out = single_run(key)

  labels = jax.nn.one_hot(out["labels"],n_clusters)
  cat = (labels * X_weight[:,None]).sum(0) / X_weight.sum()
  return {**out, "cat":cat}

def kmeans(x, x_weights, k, seed=0, max_iter=300):
  N,L,A = x.shape
  if k == 1:
    kms = {"means":(x*x_weights[:,None,None]).sum(0,keepdims=True)/x_weights.sum(),
           "labels":jnp.zeros(N,dtype=int),
           "cat":jnp.ones((1,))}
  else:
    kms = _kmeans(x.reshape(N,-1), x_weights, n_clusters=k, max_iter=max_iter, seed=seed)
    kms["means"] = kms["means"].reshape(k,L,A)
  return kms

def kmeans_sample(msa, msa_weights, k=1, samples=None, seed=0):

  assert k > 0
  
  # run kmeans
  kms = kmeans(jnp.asarray(msa), jnp.asarray(msa_weights), k=k, seed=seed)

  # sample sequences from kmeans
  key = jax.random.PRNGKey(seed)
  N,L,A = msa.shape
  if samples is None:
    # if number of samples is undefined, set to size of input MSA
    samples = N
    sampled_labels = kms["labels"]
  else:
    # sample labels
    key, key_ = jax.random.split(key)
    sampled_labels = jnp.sort(jax.random.choice(key_,jnp.arange(k),shape=(samples,),p=kms["cat"]))

  # sample MSA
  sampled_msa = kms["means"][sampled_labels]
  sampled_msa = (sampled_msa.cumsum(-1) >= jax.random.uniform(key, shape=(samples,L,1))).argmax(-1)
  o = {"kms":kms,
       "sampled_labels":sampled_labels,
       "sampled_msa":sampled_msa}
  
  return jax.tree_map(lambda x:np.asarray(x),o)