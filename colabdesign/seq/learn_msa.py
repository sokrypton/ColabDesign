import jax
import jax.numpy as jnp
import numpy

from colabdesign.seq.kmeans import kmeans
from colabdesign.seq.stats import get_stats, get_eff

# LEARN SEQUENCES
# "parameter-free" model, where we learn msa to match the statistics.
# We can take kmeans to the "next" level and directly optimize sequences to match desired stats.

class LEARN_MSA:
  def __init__(self, X, X_weight=None, samples=None,               
               mode="tied", k=1,
               seed=0, learning_rate=1e-3):
    
    assert mode in ["tied","full"]
    assert k > 0

    key = jax.random.PRNGKey(seed)
    self.k = k

    # collect X stats 
    N,L,A = X.shape
    if samples is None: samples = N
    X = jnp.asarray(X)
    X_weight = get_eff(X) if X_weight is None else jnp.asarray(X_weight)

    # run kmeans
    self.kms = kmeans(X, X_weight, k=self.k)
    stats_args = dict(add_f_ij=True, add_mf_ij=(mode=="full"), add_c=True)
    self.X_stats = get_stats(X, X_weight, labels=jax.nn.one_hot(self.kms["labels"],self.k), **stats_args)

    if samples == N:
      self.Y_labels = self.kms["labels"]
    else:
      # sample labels
      key,key_ = jax.random.split(key)
      self.Y_labels = jnp.sort(jax.random.choice(key_, jnp.arange(k), shape=(samples,), p=self.kms["cat"]))

    key, key_ = jax.random.split(key)
    Neff = X_weight.sum()
    Y_logits = jnp.log(self.kms["means"] * Neff + 0.01 * jnp.log(Neff))[self.Y_labels]
    Y = jax.nn.softmax(Y_logits + jax.random.gumbel(key_,(samples,L,A)))

    # setup the model
    def model(params, X_stats):
      # categorical reparameterization of Y
      Y_hard = jax.nn.one_hot(params["Y"].argmax(-1),A)
      Y = jax.lax.stop_gradient(Y_hard - params["Y"]) + params["Y"]

      # collect Y stats
      Y_stats = get_stats(Y, labels=jax.nn.one_hot(self.Y_labels, self.k), **stats_args)
      
      # define loss function
      i,ij = ("f_i","c_ij") if k == 1 else ("mf_i",("c_ij" if mode == "tied" else "mc_ij"))
      loss_i = jnp.square(X_stats[i] - Y_stats[i]).sum((-1,-2))
      loss_ij = jnp.square(X_stats[ij] - Y_stats[ij]).sum((-1,-2,-3)).mean(-1)
      
      if self.k > 1:
        loss_i = (loss_i * self.kms["cat"]).sum()
        if mode == "full":
          loss_ij = (loss_ij * self.kms["cat"]).sum()
      
      loss = loss_i + loss_ij

      aux = {"r":get_r(X_stats["c_ij"], Y_stats["c_ij"])}
      return loss, aux

    # setup optimizer
    self.n = 0
    init_fun, self.update_fun, self.get_params = adam(learning_rate)
    self.state = init_fun({"Y":Y})
    self.grad = jax.jit(jax.value_and_grad(model, has_aux=True))

  def get_msa(self):
    Y = np.array(self.get_params(self.state)["Y"])
    return {"kms":self.kms,
            "sampled_msa":Y.argmax(-1),
            "sampled_labels":self.Y_labels}
      
  def fit(self, steps=100, verbose=True):
    '''train model'''
    for n in range(steps):
      (loss, aux), grad = self.grad(self.get_params(self.state), self.X_stats)
      self.state = self.update_fun(self.n, grad, self.state)
      self.n += 1
      if (n+1) % (steps // 10) == 0:
        print(self.n, loss, aux["r"])