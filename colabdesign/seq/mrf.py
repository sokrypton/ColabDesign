############################
# TODO: remove reference to laxy, clean up the code
############################

def sample_msa(samples=10000, burn_in=1, temp=1.0,
               order=None, ar=False, diff=False, seq=True):

  def sample_cat(key, logits=None, probs=None):
    if logits is not None:
      hard = jax.nn.one_hot(jax.random.categorical(key,logits/temp),logits.shape[-1])
      probs = jax.nn.softmax(logits,-1)
    elif probs is not None:
      hard = (probs.cumsum(-1) >= jax.random.uniform(key, shape=probs.shape[:-1])).argmax(-1)
    if diff: hard = jax.lax.stop_gradient(hard - probs) + probs
    return hard

  def sample_pll(key,msa,par):
    N,L,A = msa.shape    
    if seq and ("w" in par or "mw" in par):
      # sequential sampling
      # burn_in = 1: autoregressive
      # burn_in > 1: gibbs
      def loop(m,x):
        i,k = x
        m_logits = []
        if "w" in par: m_logits.append(jnp.einsum("njb,ajb->na", m, par["w"][i]))
        if "b" in par: m_logits.append(par["b"][i])
        if "mw" in par: m_logits.append(jnp.einsum("nc,njb,cajb->na", par["labels"], m, par["mw"][:,i]))
        if "mb" in par: m_logits.append(jnp.einsum("nc,ca->na", par["labels"], par["mb"][:,i]))
        return m.at[:,i].set(sample_cat(k,sum(m_logits))), None
      # scan over positions
      if order is not None: i = order
      elif ar: i = jnp.arange(L)
      else: i = jax.random.permutation(key,jnp.arange(L))
      k = jax.random.split(key,L)
      return jax.lax.scan(loop,msa,(i,k))[0]
    else:
      # sample all position independently
      logits = []
      if "b" in par: logits.append(par["b"])
      if "mb" in par: logits.append(jnp.einsum("nc,cia->nia", par["labels"], par["mb"]))
      if burn_in > 1:
        if "w" in par: logits.append(jnp.einsum("njb,iajb->nia", msa, par["w"]))
        if "mw" in par: logits.append(jnp.einsum("nc,njb,ciajb->nia", par["labels"], msa, par["mw"]))
      return sample_cat(key,sum(logits))

  def sample(key, params):
    for p in ["b","w","mb","mw"]:
      if p in params:
        L,A = params[p].shape[-2:]
        break
    msa = jnp.zeros((samples,L,A))

    # sample from mixture
    if "c" in params and ("mb" in params or "mw" in params):
      c = params["c"]
      labels_logits = jnp.tile(c,(samples,1))
      params["labels"] = jax.nn.one_hot(jax.random.categorical(key,labels_logits),c.shape[0])
      if diff:
        labels_soft = jax.nn.softmax(labels_logits,-1)
        params["labels"] = jax.lax.stop_gradient(params["labels"] - labels_soft) + labels_soft
    else:
      params["labels"] = None

    # number of iterations (burn-in)
    sample_loop = lambda m,k:(sample_pll(k,m,params),None)
    iters = jax.random.split(key,burn_in)
    msa = jax.lax.scan(sample_loop,msa,iters)[0]
    return {"msa":msa, "labels":params["labels"]}

  return sample

def reg_loss(params, lam):
  reg_loss = []
  if "b" in params:
    reg_loss.append(lam * jnp.square(params["b"]).sum())
  if "w" in params:
    L,A = params["w"].shape[-2:]
    reg_loss.append(lam/2*(L-1)*(A-1) * jnp.square(params["w"]).sum())
  if "mb" in params:
    reg_loss.append(lam * jnp.square(params["mb"]).sum())
  if "mw" in params:
    L,A = params["mw"].shape[-2:]
    reg_loss.append(lam/2*(L-1)*(A-1) * jnp.square(params["mw"]).sum())
  return sum(reg_loss)

def pll_loss(params, inputs, order=None, labels=None):
  logits = []

  L = inputs["x"].shape[1]
  w_mask = 1-jnp.eye(L)
  if order is not None:
    w_mask *= ar_mask(order)

  if "b" in params:
    logits.append(params["b"])
  if "w" in params:
    w = params["w"]
    w = 0.5 * (w + w.transpose([2,3,0,1])) * w_mask[:,None,:,None]
    logits.append(jnp.einsum("nia,iajb->njb", inputs["x"], w))

  # MIXTURES
  if "mb" in params:
    logits.append(jnp.einsum("nc,cia->nia", labels, params["mb"]))

  if "mw" in params:
    mw = params["mw"]
    mw = 0.5 * (mw + mw.transpose([0,3,4,1,2])) * w_mask[None,:,None,:,None]
    logits.append(jnp.einsum("nc,nia,ciajb->njb", labels, inputs["x"], mw))
      
  # categorical-crossentropy (or pseudo-likelihood)
  cce_loss = -(inputs["x"] * jax.nn.log_softmax(sum(logits))).sum([1,2])

  return (cce_loss*inputs["x_weight"]).sum()

class MRF:  
  def __init__(self, X, X_weight=None,
               batch_size=None,
               ar=False, ar_ent=False,
               lam=0.01,
               k=1, lr=0.1, shared=False, tied=True, full=False):
    
    ## MODE ##
    inc = ["b","w"] if (tied or full) else ["b"]
    if k > 1:
      if shared:
        if tied: inc += ["mb"]
        if full: inc += ["mb","mw"]
      else:
        if tied: inc = ["mb","w"]
        if full: inc = ["mb","mw"]

    N,L,A = X.shape
    self.batch_size = batch_size
    self.k = k

    # weight per sequence
    X = jnp.asarray(X)
    X_weight = get_eff(X) if X_weight is None else jnp.asarray(X_weight)
    self.Neff = X_weight.sum()

    if batch_size is None:
      learning_rate = lr * np.log(N)/L
    else:
      lam = lam * batch_size/N
      learning_rate = lr * jnp.log(batch_size)/L

    if ar:
      self.order = jnp.arange(L)
    elif ar_ent:
      f_i = (X * X_weight[:,None,None]).sum(0)/self.Neff
      self.order = (-f_i * jnp.log(f_i + 1e-8)).sum(-1).argsort()
    else:
      self.order = None

    # setup the model
    def model(params, inputs):
      labels = inputs["labels"] if "labels" in inputs else None
      pll = pll_loss(params, inputs, self.order, labels)
      reg = reg_loss(params, lam)
      loss = pll + reg
      return None, loss

    # initialize inputs
    self.inputs = {"x":X, "x_weight":X_weight}

    # initialize params
    self.params = {}
    if "w" in inc: self.params["w"] = jnp.zeros((L,A,L,A))
    if "mw" in inc: self.params["mw"] = jnp.zeros((k,L,A,L,A))
    if "b" in inc:
      b = jnp.log((X * X_weight[:,None,None]).sum(0) + (lam+1e-8) * jnp.log(self.Neff))
      self.params["b"] = b - b.mean(-1,keepdims=True)
    if "mb" in inc or "mw" in inc:
      kms = kmeans(X, X_weight, k=k)
      self.inputs["labels"] = kms["labels"]
      mb = jnp.log(kms["means"] * self.Neff + (lam+1e-8) * jnp.log(self.Neff))
      self.params["mb"] = mb - mb.mean(-1,keepdims=True)
      if "b" in self.params:
        self.params["mb"] -= self.params["b"]

    # setup optimizer
    self.opt = laxy.OPT(model, self.params, lr=learning_rate)

  def get_msa(self, samples=1000, burn_in=1):
    self.params = self.opt.get_params()
    if "labels" in self.inputs:
      self.params["c"] = jnp.log((self.inputs["x_weight"][:,None] * self.inputs["labels"]).sum(0) + 1e-8)
    
    key = laxy.get_random_key()
    return sample_msa(samples=samples,burn_in=burn_in,order=self.order)(key, self.params)
  
  def get_w(self):
    self.params = self.opt.get_params()
    w = []
    if "w" in self.params: w.append(self.params["w"])
    if "mw" in self.params: w.append(self.params["mw"].sum(0))
    w = sum(w)
    w = (w + w.transpose(2,3,0,1))/2
    w = w - w.mean((1,3),keepdims=True)
    return w
    
  def fit(self, steps=100, verbose=True, return_losses=False):
    '''train model'''
    losses = self.opt.fit(self.inputs, steps=steps, batch_size=self.batch_size,
                         verbose=verbose, return_losses=return_losses)
    if return_losses: return losses

class MRF_BM:
  def __init__(self, X, X_weight=None, samples=1000,
               burn_in=1, temp=1.0,
               ar=False, ar_ent=True,
               lr=0.05, lam=0.01,
               k=1, mode="tied"):

    ## MODE ##
    inc = ["b","mb"] if k > 1 else ["b"]
    if mode == "tied": inc += ["w"]
    if mode == "full": inc += ["w","mw"] if k > 1 else ["w"]

    self.X = jnp.asarray(X)
    N,L,A = self.X.shape
    learning_rate = lr * np.log(N)/L

    # weight per sequence
    self.X_weight = get_eff(X) if X_weight is None else jnp.asarray(X_weight)
    self.Neff = self.X_weight.sum()

    # collect stats 
    if k > 1:
      self.kms = kmeans(self.X, self.X_weight, k=k)
      self.labels = self.kms["labels"]
      self.inputs = get_stats(self.X, self.X_weight, labels=self.kms["labels"],
                              add_mf_ij=("mw" in inc))
      self.inputs["c"] = self.kms["cat"]
    else:
      self.labels = None
      self.inputs = get_stats(self.X, self.X_weight)
    
    # low entropy to high entropy
    if ar_ent:
      ent = -(self.inputs["f_i"] * jnp.log(self.inputs["f_i"] + 1e-8)).sum(-1)
      self.order = ent.argsort()
    elif ar: self.order = jnp.arange(L)
    else: self.order = None

    self.burn_in = burn_in
    self.temp = temp

    # setup the model
    def model(params, inputs):

      # sample msa
      sample = sample_msa(samples=samples, burn_in=burn_in,
                          temp=temp, order=self.order)(inputs["key"], params)

      # compute stats
      stats = get_stats(sample["msa"], labels=sample["labels"],
                        add_mf_ij=("mw" in params))

      # define gradients
      grad = {}
      I = (1-jnp.eye(L))[:,None,:,None]
      if "c" in params: grad["c"] = sample["labels"].mean(0) - inputs["c"]
      if "b" in params: grad["b"] = stats["f_i"] - inputs["f_i"]
      if "w" in params: grad["w"] = (stats["f_ij"] - inputs["f_ij"]) * I
      if "mb" in params: grad["mb"] = stats["mf_i"] - inputs["mf_i"]
      if "mw" in params: grad["mw"] = (stats["mf_ij"] - inputs["mf_ij"]) * I[None]

      # add regularization
      reg_grad = jax.grad(reg_loss)(params,lam)

      for g in grad.keys():
        if g in reg_grad: grad[g] = grad[g] * self.Neff + reg_grad[g]
        else: grad[g] = grad[g] * self.Neff

      return None, None, grad

    # initialize model params
    self.params = {}
    if "w" in inc:
      self.params["w"] = jnp.zeros((L,A,L,A))
    if "b" in inc:
      b = jnp.log(self.inputs["f_i"] * self.Neff + (lam+1e-8) * jnp.log(self.Neff))
      self.params["b"] = b - b.mean(-1,keepdims=True)

    # setup mixture params
    if "mb" in inc or "mw" in inc:
      c = jnp.log(self.inputs["c"] * self.Neff + 1e-8)
      self.params["c"] = c - c.mean(-1,keepdims=True)

    if "mw" in inc:
      self.params["mw"] = jnp.zeros((k,L,A,L,A))
    if "mb" in inc:
      mb = jnp.log(self.inputs["mf_i"] * self.Neff + (lam+1e-8) * jnp.log(self.Neff))
      self.params["mb"] = mb - mb.mean(-1,keepdims=True)
      if "b" in self.params: self.params["mb"] -= self.params["b"]

    # setup optimizer
    self.opt = laxy.OPT(model, self.params, lr=learning_rate, has_grad=True)
  
  def get_msa(self, samples=1000, burn_in=None, temp=None, seed=0):
    if burn_in is None: burn_in = self.burn_in
    if temp is None: temp = self.temp

    self.params = self.opt.get_params()
    key = jax.random.PRNGKey(seed)
    return sample_msa(samples=samples,
                      burn_in=self.burn_in,
                      order=self.order)(key, self.params)["msa"]

  def fit(self, steps=1000, verbose=True):
    '''train model'''
    self.opt.fit(self.inputs, steps=steps, verbose=verbose)