import random, copy, os
import jax
import jax.numpy as jnp
import numpy as np
from alphafold.common import residue_constants
from colabdesign.af.utils import update_dict

try:
  from jax.example_libraries.optimizers import sgd, adam
except:
  from jax.experimental.optimizers import sgd, adam

####################################################
# AF_DESIGN - design functions
####################################################
class _af_design:
  def _setup_optimizer(self, optimizer="sgd", lr_scale=1.0, **kwargs):
    '''setup which optimizer to use'''
    if optimizer == "adam":
      optimizer = adam
      lr = 0.02 * lr_scale
    elif optimizer == "sgd":
      optimizer = sgd
      lr = 0.1 * lr_scale
    else:
      lr = lr_scale
      
    self._init_fun, self._update_fun, self._get_params = optimizer(lr)
    self._k = 0

  def _init_seq(self, mode=None, seq=None, rm_aa=None, add_seq=False, **kwargs):
    '''initialize sequence'''

    # backward compatibility
    if "seq_init" in kwargs:
      x = kwargs["seq_init"]
      if isinstance(x, np.ndarray) or isinstance(x, jnp.ndarray): seq = x
      elif isinstance(x, str):
        # TODO this could result in some issues...
        if len(x) == self._len: seq = x
        else: mode = x

    if mode is None: mode = ""
    
    # initialize sequence
    shape = (self.args["num_seq"], self._len, 20)

    self._key, key = jax.random.split(self._key)
    x = 0.01 * jax.random.normal(key, shape)

    # initialize bias
    bias = jnp.zeros(shape[1:])
    use_bias = False

    if "gumbel" in mode:
      self._key, key = jax.random.split(self._key)
      x = jax.random.gumbel(key, shape) / 2.0

    if "soft" in mode:
      x = jax.nn.softmax(x * 2.0)

    if "wildtype" in mode or "wt" in mode:
      wt = jax.nn.one_hot(self._wt_aatype, 20)
      if self.protocol == "fixbb":
        wt = wt[...,:self._len]
      if self.protocol == "partial":
        x = x.at[...,self.opt["pos"],:].set(wt)
        if add_seq:
          bias = bias.at[self.opt["pos"],:].set(wt * 1e8)
          use_bias = True
      else:
        x = x.at[...,:,:].set(wt)
        if add_seq:
          bias = bias.at[:,:].set(wt * 1e8)
          use_bias = True

    if seq is not None:
      if isinstance(seq, np.ndarray) or isinstance(seq, jnp.ndarray):
        x_ = seq
      elif isinstance(seq, str):
        x_ = jnp.array([residue_constants.restype_order.get(aa,-1) for aa in seq])
        x_ = jax.nn.one_hot(x_,20)
      x = x + x_
      if add_seq:
        bias = bias + x_ * 1e8
        use_bias = True 

    if rm_aa is not None:
      for aa in rm_aa.split(","):
        bias = bias - 1e8 * np.eye(20)[residue_constants.restype_order[aa]]
        use_bias = True

    if use_bias:
      self.opt["bias"] = bias
    self.params = {"seq":x}
    self._state = self._init_fun(self.params)


  def set_opt(self, *args, **kwargs):
    update_dict(self.opt, *args, **kwargs)
    
  def set_weights(self, *args, **kwargs):
    update_dict(self.opt["weights"], *args, **kwargs)

  def restart(self, seed=None, weights=None, opt=None, set_defaults=False, keep_history=False, **kwargs):
    
    # set weights and options
    if set_defaults:
      update_dict(self._default_opt, opt)
      update_dict(self._default_opt["weights"], weights)
        
      if not self.use_struct:
        # remove structure module specific weights
        struct_list = ["rmsd","fape","plddt","pae"]
        keys = list(self._default_opt["weights"].keys())
        for k in keys:
          if k.split("_")[-1] in struct_list:
            self._default_opt["weights"].pop(k)
    
    if not keep_history:
      self.opt = copy.deepcopy(self._default_opt)

    self.set_opt(opt)
    self.set_weights(weights)

    # setup optimizer
    self._setup_optimizer(**kwargs)
    
    # initialize sequence
    self._seed = random.randint(0,2147483647) if seed is None else seed
    self._key = jax.random.PRNGKey(self._seed)
    self._init_seq(**kwargs)

    # initialize trajectory
    if not keep_history:
      if self.use_struct:
        self._traj = {"losses":[],"seq":[],"xyz":[],"plddt":[],"pae":[]}
      else:
        self._traj = {"losses":[],"seq":[],"con":[]}
      self._best_loss, self._best_aux = np.inf, None
      self._best_outs = self._best_aux
  
  def clear_best(self):
    self._best_loss, self._best_aux = p.inf, None
    self._best_outs = self._best_aux
  
  def _save_results(self, save_best=False, verbose=True):
    '''save the results and update trajectory'''

    # save best result
    if save_best and self.loss < self._best_loss:
      self._best_loss = self.loss
      self._best_aux = self.aux
      # backward compatibility
      self._best_outs = self._best_aux

    # save losses
    losses = jax.tree_map(float, self.aux["losses"])
    losses.update({"models":self.aux["model_num"],
                   "recycles":int(self.aux["recycles"]),
                   "loss":float(self.loss),
                   "hard":int(self.opt["hard"]),
                   "soft":float(self.opt["soft"]),
                   "temp":float(self.opt["temp"])})
    
    if self.protocol == "fixbb" or (self.protocol == "binder" and self._redesign):
      # compute sequence recovery
      _aatype = self.aux["seq"]["hard"].argmax(-1)
      L = min(_aatype.shape[-1], self._wt_aatype.shape[-1])
      losses["seqid"] = float((_aatype[...,:L] == self._wt_aatype[...,:L]).mean())

    # print losses  
    if verbose:
      SEEN = []
      
      def to_str(ks):
        out = []
        for k in ks:
          SEEN.append(k)
          if k in losses and ("rmsd" in k or self.opt["weights"].get(k,True)):
            out.append(f"{k}")
            v = losses[k]
            if isinstance(v,int): f = 0
            elif isinstance(v,float): f = 2
            else: f = None
            out.append(f"{v}" if f is None else f"{v:.{f}f}")
        return out
      
      out = [to_str(["models","recycles",
                     "hard","soft","temp","seqid","loss",
                     "msa_ent","plddt","pae","helix","con",
                     "i_pae","i_con",
                     "sc_fape","sc_rmsd",
                     "dgram_cce","fape","6D","rmsd"])]

      for k,v in losses.items():
        if k not in SEEN: out.append(to_str([k]))

      print_str = " ".join(sum(out,[]))
      print(f"{self._k}\t{print_str}")

    # save trajectory
    traj = {"losses":losses, "seq":np.asarray(self.aux["seq"]["pseudo"])}
    if self.use_struct:
      traj["xyz"] = np.asarray(self.aux["final_atom_positions"][:,1,:])
      traj["plddt"] = np.asarray(self.aux["plddt"])
      if "pae" in self.aux: traj["pae"] = np.asarray(self.aux["pae"])
    else:
      traj["con"] = np.asarray(self.aux["contact_map"])
    for k,v in traj.items(): self._traj[k].append(v)
    
  def _recycle(self, model_params, backprop=True):                
    # initialize previous (for recycle)
    L = self._inputs["residue_index"].shape[-1]
    
    if self.use_struct:
      self._inputs["prev"] = {'prev_msa_first_row': np.zeros([L,256]),
                              'prev_pair': np.zeros([L,L,128]),
                              'prev_pos': np.zeros([L,37,3])}
    else:
      self._inputs["prev"] = {'prev_msa_first_row': np.zeros([L,256]),
                              'prev_pair': np.zeros([L,L,128]),
                              'prev_dgram': np.zeros([L,L,64])}
              
    if self.args["recycle_mode"] in ["backprop","add_prev"]:
      self.opt["recycles"] = self._runner.config.model.num_recycle

    recycles = self.opt["recycles"]

    if self.args["recycle_mode"] == "average":
      # average gradient across recycles
      if backprop:
        grad = []
      for r in range(recycles+1):
        self._key, key = jax.random.split(self._key)
        if backprop:
          (loss, aux), _grad = self._grad_fn(self.params, model_params, self._inputs, key, self.opt)
          grad.append(_grad)
        else:
          loss, aux = self._fn(self.params, model_params, self._inputs, key, self.opt)
        self._inputs["prev"] = aux["prev"]
      if backprop:
        grad = jax.tree_map(lambda *x: jnp.asarray(x).mean(0), *grad)
    else:
      if self.args["recycle_mode"] == "sample":
        self._key, key = jax.random.split(self._key)
        self.opt["recycles"] = jax.random.randint(key, [], 0, recycles+1)
      
      self._key, key = jax.random.split(self._key)
      if backprop:
        (loss, aux), grad = self._grad_fn(self.params, model_params, self._inputs, key, self.opt)
      else:
        loss, aux = self._fn(self.params, model_params, self._inputs, key, self.opt)

    aux["recycles"] = self.opt["recycles"]
    self.opt["recycles"] = recycles
    
    if backprop:
      return (loss, aux), grad
    else:
      return loss, aux
    
  def run(self, backprop=True):
    '''run model to get outputs, losses and gradients'''
    
    # decide which model params to use
    m = self.opt["models"]
    ns = jnp.arange(2) if self.args["use_templates"] else jnp.arange(5)
    if self.args["model_sample"] and m != len(ns):
      self._key, key = jax.random.split(self._key)
      model_num = jax.random.choice(key,ns,(m,),replace=False)
    else:
      model_num = ns[:m]
    model_num = np.array(model_num).tolist()    
    
    # loop through model params
    _loss, _aux, _grad = [],[],[]
    for n in model_num:
      p = self._model_params[n]
      if backprop:
        (l,a),g = self._recycle(p)
        _grad.append(g)
      else:
        l,a = self._recycle(p, backprop=False)
      _loss.append(l); _aux.append(a)

    # update gradients
    if backprop:
      self.grad = jax.tree_map(lambda *v: jnp.stack(v).mean(0), *_grad)
    else:
      self.grad = jax.tree_map(jnp.zeros_like, self.params)

    # update loss (take mean across models)
    self.loss = jnp.stack(_loss).mean()
    
    # update aux
    _aux = jax.tree_map(lambda *v: jnp.stack(v), *_aux)    
    self.aux = jax.tree_map(lambda x:x[0], _aux) # pick one example
    self.aux["losses"] = jax.tree_map(lambda v: v.mean(0), _aux["losses"])
    self.aux["model_num"] = model_num    
    self._outs = self.aux # backward compatibility
    
  #-------------------------------------
  # STEP FUNCTION
  #-------------------------------------
  def step(self, lr_scale=1.0, backprop=True, callback=None, **kwargs):
    '''do one step of gradient descent'''
    
    # run
    self.run(backprop=backprop)
    if callback is not None: callback(self)

    # normalize gradient
    g = self.grad["seq"]
    gn = jnp.linalg.norm(g,axis=(-1,-2),keepdims=True)
    _len = self._len # TODO (we might need to change this for partial design)
    self.grad["seq"] *= lr_scale * jnp.sqrt(_len)/(gn+1e-7)

    # apply gradient
    self._state = self._update_fun(self._k, self.grad, self._state)
    self.params = self._get_params(self._state)

    self._k += 1
    self._save_results(**kwargs)


  #-----------------------------------------------------------------------------
  # DESIGN FUNCTIONS
  #-----------------------------------------------------------------------------
  def design(self, iters,
             temp=1.0, e_temp=None,
             soft=False, e_soft=None,
             hard=False, dropout=True, gumbel=False, 
             opt=None, weights=None, callback=None, **kwargs):
      
    # override settings if defined
    self.set_opt(opt)
    self.set_weights(weights)
      
    self.set_opt(hard=hard, dropout=dropout, gumbel=gumbel)
    if e_soft is None: e_soft = soft
    if e_temp is None: e_temp = temp
    for i in range(iters):
      self.opt["temp"] = e_temp + (temp - e_temp) * (1-i/(iters-1)) ** 2
      self.opt["soft"] = soft + (e_soft - soft) * i/(iters-1)
      # decay learning rate based on temperature
      lr_scale = (1 - self.opt["soft"]) + (self.opt["soft"] * self.opt["temp"])
      self.step(lr_scale=lr_scale, callback=callback, **kwargs)

  def design_logits(self, iters, **kwargs):
    '''optimize logits'''
    self.design(iters, **kwargs)

  def design_soft(self, iters, **kwargs):
    ''' optimize softmax(logits/temp)'''
    self.design(iters, soft=True, **kwargs)
  
  def design_hard(self, iters, **kwargs):
    ''' optimize argmax(logits)'''
    self.design(iters, soft=True, hard=True, **kwargs)

  def design_2stage(self, soft_iters=100, temp_iters=100, hard_iters=50,
                    temp=1.0, dropout=True, gumbel=False, **kwargs):
    '''two stage design (soft→hard)'''
    self.design(soft_iters, soft=True, temp=temp, dropout=dropout, gumbel=gumbel, **kwargs)
    self.design(temp_iters, soft=True, temp=temp, dropout=dropout, gumbel=False,  e_temp=1e-2, **kwargs)
    self.design(hard_iters, soft=True, temp=1e-2, dropout=False,   gumbel=False,  hard=True, save_best=True, **kwargs)

  def design_3stage(self, soft_iters=300, temp_iters=100, hard_iters=50, 
                    temp=1.0, dropout=True, gumbel=False, **kwargs):
    '''three stage design (logits→soft→hard)'''
    self.design(soft_iters, e_soft=True, temp=temp, dropout=dropout, gumbel=gumbel, **kwargs)
    self.design(temp_iters, soft=True,   temp=temp, dropout=dropout, gumbel=False, e_temp=1e-2,**kwargs)
    self.design(hard_iters, soft=True,   temp=1e-2, dropout=False,   gumbel=False, hard=True, save_best=True, **kwargs)  
    
  def template_predesign(self, iters=100, soft=True, hard=False, dropout=True, temp=1.0, **kwargs):
    '''use template for predesign stage, then increase template dropout until gone'''
    self.set_opt(hard=hard, soft=soft, dropout=dropout, temp=temp)
    for i in range(iters):
      self.opt["template_dropout"] = i/(iters-1)
      self.step(**kwargs)
  
  def design_semigreedy(self, iters=100, tries=20, dropout=False,
                        use_plddt=True, save_best=True, verbose=True):
    '''semigreedy search'''
    
    self.set_opt(hard=True, soft=True, temp=1.0, dropout=dropout)
    if self._k == 0:
      self.run(backprop=False)

    def mut(seq, plddt=None):
      '''mutate random position'''
      L,A = seq.shape[-2:]
      while True:
        self._key, key_L, key_A = jax.random.split(self._key,3)

        if plddt is None:
          i = jax.random.randint(key_L,[],0,L)
        else:
          p = (1-plddt)/(1-plddt).sum(-1,keepdims=True)
          # bias mutations towards positions with low pLDDT
          # https://www.biorxiv.org/content/10.1101/2021.08.24.457549v1
          i = jax.random.choice(key_L,jnp.arange(L),[],p=p)

        a = jax.random.randint(key_A,[],0,A)
        if seq[0,i,a] == 0: break      
      return seq.at[:,i,:].set(jnp.eye(A)[a])

    def get_seq():
      return jax.nn.one_hot(self.params["seq"].argmax(-1),20)

    # optimize!
    for i in range(iters):
      seq = get_seq()
      plddt = None
      if use_plddt:
        plddt = np.asarray(self.aux["plddt"])
        if self.protocol == "binder":
          plddt = plddt[self._target_len:]
        else:
          plddt = plddt[:self._len]
      
      buff = []
      for _ in range(tries):
        self.params["seq"] = mut(seq, plddt)
        self.run(backprop=False)
        buff.append({"aux":self.aux, "loss":self.loss, "seq":self.params["seq"]})
      
      # accept best      
      buff = buff[np.argmin([x["loss"] for x in buff])]
      self.aux, self.loss, self.params["seq"] = buff["aux"], buff["loss"], buff["seq"]
      self._k += 1
      self._save_results(save_best=save_best, verbose=verbose)