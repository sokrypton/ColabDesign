import random, copy, os
import jax
import jax.numpy as jnp
import numpy as np

from alphafold.common import residue_constants

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
    self._key, _key = jax.random.split(self._key)
    shape = (self.args["num_seq"], self._len, 20)
    x = 0.01 * jax.random.normal(_key, shape)

    # initialize bias
    b = jnp.zeros(shape[1:])

    if "gumbel" in mode:
      x = jax.random.gumbel(_key, shape) / 2.0

    if "soft" in mode:
      x = jax.nn.softmax(x * 2.0)

    if "wildtype" in mode or "wt" in mode:
      wt = jax.nn.one_hot(self._wt_aatype, 20)
      if self.protocol == "fixbb":
        wt = wt[...,:self._len]
      elif self.protocol == "partial":
        x = x.at[...,self.opt["pos"],:].set(wt)
        if add_seq: b = b.at[self.opt["pos"],:].set(wt * 1e8)
      else:
        x = x.at[...,:,:].set(wt)
        if add_seq: b = b.at[:,:].set(wt * 1e8)

    if seq is not None:
      if isinstance(seq, np.ndarray) or isinstance(seq, jnp.ndarray):
        x_ = seq
      elif isinstance(seq, str):
        x_ = jnp.array([residue_constants.restype_order.get(aa,-1) for aa in seq])
        x_ = jax.nn.one_hot(x_,20)
      x = x + x_
      if add_seq: b = b + x_ * 1e8

    if rm_aa is not None:
      for aa in rm_aa.split(","):
        b = b - 1e8 * np.eye(20)[residue_constants.restype_order[aa]]

    self.opt["bias"] = b
    self._params = {"seq":x}
    self._state = self._init_fun(self._params)

  def restart(self, seed=None, weights=None, opt=None, set_defaults=False, keep_history=False, **kwargs):
    
    # set weights and options
    if set_defaults:
      if weights is not None: self._default_weights.update(weights)
      if opt is not None: self._default_opt.update(opt)
        
      if not self.use_struct:
        # remove structure module specific weights
        struct_list = ["rmsd","fape","plddt","pae"]
        keys = list(self._default_weights.keys())
        for k in keys:
          if k.split("_")[-1] in struct_list:
            self._default_weights.pop(k)
    
    if not keep_history:
      self.opt = copy.deepcopy(self._default_opt)
      self.opt["weights"] = self._default_weights.copy()

    if weights is not None: self.opt["weights"].update(weights)
    if opt is not None: self.opt.update(opt)

    # setup optimizer
    self._setup_optimizer(**kwargs)
    
    # initialize sequence
    self._seed = random.randint(0,2147483647) if seed is None else seed
    self._key = jax.random.PRNGKey(self._seed)
    self._init_seq(**kwargs)

    # initialize trajectory
    if not keep_history:
      if self.use_struct:
        self._traj = {"seq":[],"xyz":[],"plddt":[],"pae":[]}
      else:
        self._traj = {"seq":[],"con":[]}
      self.losses, self._best_loss, self._best_outs = [], np.inf, None
      
  def _save_results(self, save_best=False, verbose=True):
    '''save the results and update trajectory'''

    # save best result
    if save_best and self._loss < self._best_loss:
      self._best_loss = self._loss
      self._best_outs = self._outs
    
    # compile losses
    self._losses.update({"model":self._model_num,"loss":self._loss, 
                        **{k:self.opt[k] for k in ["soft","hard","temp"]}})
    
    if self.protocol == "fixbb" or (self.protocol == "binder" and self._redesign):
      # compute sequence recovery
      _aatype = self._outs["seq"]["hard"].argmax(-1)
      L = min(_aatype.shape[-1], self._wt_aatype.shape[-1])
      self._losses["seqid"] = (_aatype[...,:L] == self._wt_aatype[...,:L]).mean()

    # save losses
    self.losses.append(self._losses)

    # print losses  
    if verbose:
      def to_str(ks, f=2):
        out = []
        for k in ks:
          if k in self._losses and ("rmsd" in k or self.opt["weights"].get(k,True)):
            out.append(f"{k}")
            if f is None: out.append(f"{self._losses[k]}")
            else: out.append(f"{self._losses[k]:.{f}f}")
        return out
      out = [to_str(["model","recycles"],None),
             to_str(["soft","temp","seqid","loss"]),
             to_str(["msa_ent","plddt","pae","helix","con",
                     "i_pae","i_con",
                     "sc_fape","sc_rmsd",
                     "dgram_cce","fape","6D","rmsd"])]
      print_str = " ".join(sum(out,[]))
      print(f"{self._k}\t{print_str}")

    # save trajectory
    traj = {"seq":self._outs["seq"]["pseudo"]}
    if self.use_struct:
      traj.update({"xyz":self._outs["final_atom_positions"][:,1,:],
                   "plddt":self._outs["plddt"]})
      if "pae" in self._outs:
        traj.update({"pae":self._outs["pae"]})
    else:
      traj["con"] = self._outs["contact_map"]
      
    for k,v in traj.items(): self._traj[k].append(np.array(v))
    
  def _recycle(self, model_params, key, params=None, opt=None):
    if params is None: params = self._params
    if opt is None: opt = self.opt
    if self.args["recycle_mode"] == "average":
      #----------------------------------------
      # average gradients across all recycles
      #----------------------------------------
      L = self._inputs["residue_index"].shape[-1]
      self._inputs.update({'init_pos': np.zeros([1, L, 37, 3]),
                           'init_msa_first_row': np.zeros([1, L, 256]),
                           'init_pair': np.zeros([1, L, L, 128])})
      
      grad = []
      for r in range(opt["recycles"]+1):
        key, _key = jax.random.split(key)
        _opt = copy.deepcopy(opt)
        for k,x in _opt["weights"].items():
          if isinstance(x, np.ndarray) or isinstance(x, jnp.ndarray) or isinstance(x, list):
            # use different weights at each recycle
            _opt["weights"][k] = float(x[r])
          else:
            _opt["weights"][k] = float(x)
        (loss,outs),_grad = self._grad_fn(params, model_params, self._inputs, _key, _opt)
        grad.append(_grad)
        self._inputs.update(outs["init"])
      grad = jax.tree_util.tree_map(lambda *x: jnp.asarray(x).mean(0), *grad)
      return (loss, outs), grad
    else:
      #----------------------------------------
      # use gradients from last recycle
      #----------------------------------------
      if self.args["recycle_mode"] == "sample":
        # decide number of recycles to use
        key, _key = jax.random.split(key)
        opt["recycles"] = int(jax.random.randint(_key,[],0,self._default_opt["recycles"]+1))
      
      return self._grad_fn(params, model_params, self._inputs, key, opt)
      
  def _update_grad(self):
    '''update the gradient'''

    # get current params
    self._params = self._get_params(self._state)

    # update key
    self._key, key = jax.random.split(self._key)

    # decide which model params to use
    m = self.opt["models"]
    ns = jnp.arange(2) if self.args["use_templates"] else jnp.arange(5)
    if self.args["model_mode"] == "fixed" or m == len(ns):
      self._model_num = ns[:m]
    elif self.args["model_mode"] == "sample":
      key,_key = jax.random.split(key)
      self._model_num = jax.random.choice(_key,ns,(m,),replace=False)
    
    #------------------------
    # run in parallel
    #------------------------
    if self.args["model_parallel"]:
      p = self._model_params
      if self.args["model_mode"] == "sample":
        p = self._sample_params(p, self._model_num)
      (l,o),g = self._recycle(p, key)
      self._grad = jax.tree_map(lambda x: x.mean(0), g)
      self._loss, self._outs = l.mean(),jax.tree_map(lambda x:x[0],o)
      self._losses = jax.tree_map(lambda x: x.mean(0), o["losses"])

    #------------------------
    # run in serial
    #------------------------
    else:
      self._model_num = np.array(self._model_num).tolist()    
      _loss, _losses, _outs, _grad = [],[],[],[]
      for n in self._model_num:
        p = self._model_params[n]
        (l,o),g = self._recycle(p, key)
        _loss.append(l); _outs.append(o); _grad.append(g)
        _losses.append(o["losses"])
      self._grad = jax.tree_util.tree_map(lambda *v: jnp.asarray(v).mean(0), *_grad)      
      self._loss, self._outs = jnp.mean(jnp.asarray(_loss)), _outs[0]
      self._losses = jax.tree_util.tree_map(lambda *v: jnp.asarray(v).mean(), *_losses)

  #-------------------------------------
  # STEP FUNCTION
  #-------------------------------------
  def _step(self, weights=None, lr_scale=1.0, **kwargs):
    '''do one step'''

    if weights is not None: self.opt["weights"].update(weights)
    
    # update gradient
    self._update_grad()

    # normalize gradient
    g = self._grad["seq"]
    gn = jnp.linalg.norm(g,axis=(-1,-2),keepdims=True)
    self._grad["seq"] *= lr_scale * jnp.sqrt(self._len)/(gn+1e-7)

    # apply gradient
    self._state = self._update_fun(self._k, self._grad, self._state)
    self._k += 1
    self._save_results(**kwargs)


  #-----------------------------------------------------------------------------
  # DESIGN FUNCTIONS
  #-----------------------------------------------------------------------------
  def design(self, iters,
              temp=1.0, e_temp=None,
              soft=False, e_soft=None,
              hard=False, dropout=True, gumbel=False, 
              opt=None, weights=None, **kwargs):

    if opt is not None: self.opt.update(opt)
    if weights is not None: self.opt["weights"].update(weights)
      
    self.opt.update({"hard":float(hard),"dropout":dropout,"gumbel":gumbel})
    if e_soft is None: e_soft = soft
    if e_temp is None: e_temp = temp
    for i in range(iters):
      self.opt["temp"] = e_temp + (temp - e_temp) * (1-i/(iters-1)) ** 2
      self.opt["soft"] = soft + (e_soft - soft) * i/(iters-1)
      # decay learning rate based on temperature
      lr_scale = (1 - self.opt["soft"]) + (self.opt["soft"] * self.opt["temp"])
      self._step(lr_scale=lr_scale, **kwargs)

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
    self.opt.update({"hard":float(hard), "soft":float(soft),
                     "dropout":dropout, "temp":temp})
    for i in range(iters):
      self.opt["template_dropout"] = i/(iters-1)
      self._step(**kwargs)
