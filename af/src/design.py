import random, copy, os
import jax
import jax.numpy as jnp
import numpy as np

from alphafold.common import residue_constants

# borrowing update_dict function from TrDesign
from tr.src.utils import update_dict

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

  def set_weights(self, *args, **kwargs):
    update_dict(self.opt["weights"], *args, **kwargs)

  def set_opt(self, *args, **kwargs):
    update_dict(self.opt, *args, **kwargs)

  def restart(self, seed=None, weights=None, opt=None, set_defaults=False, keep_history=False, **kwargs):
    
    # set weights and options
    if set_defaults:
      update_dict(self._default_weights, weights)
      update_dict(self._default_opt, opt)
        
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

    self.set_weights(weights)
    self.set_opt(opt)

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
    self._losses.update({"models":self._model_num,
                         "recycles":self._outs["recycles"],
                         "loss":self._loss,
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
      SEEN = []
      def to_str(ks, f=2):
        out = []
        for k in ks:
          SEEN.append(k)
          if k in self._losses and ("rmsd" in k or self.opt["weights"].get(k,True)):
            out.append(f"{k}")
            if f is None: out.append(f"{self._losses[k]}")
            else: out.append(f"{self._losses[k]:.{f}f}")
        return out
      
      out = [to_str(["models","recycles"],None),
             to_str(["hard"],0),
             to_str(["soft","temp","seqid","loss"]),
             to_str(["msa_ent","plddt","pae","helix","con",
                     "i_pae","i_con",
                     "sc_fape","sc_rmsd",
                     "dgram_cce","fape","6D","rmsd"])]

      for k in self._losses.keys():
        if k not in SEEN: out.append(to_str([k]))

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
    
  def run(self, seq=None, params=None, opt=None, weights=None,
          backprop=True, callback=None):
    '''run model to get outputs, losses and gradients'''
    
    # override settings if defined
    update_dict(self.params, seq=seq)
    update_dict(self.params, params)
    self.set_opt(opt)
    self.set_weights(weights)

    # decide which model params to use
    m = self.opt["models"]
    ns = jnp.arange(2) if self.args["use_templates"] else jnp.arange(5)
    if self.args["model_sample"] and m != len(ns):
      self._key, key = jax.random.split(self._key)
      self._model_num = jax.random.choice(key,ns,(m,),replace=False)
    else:
      self._model_num = ns[:m]
    
    self._model_num = np.array(self._model_num).tolist()    
    _loss, _losses, _outs, _grad = [],[],[],[]
    for n in self._model_num:
      p = self._model_params[n]
      if backprop:
        (l,o),g = self._recycle(p)
        _grad.append(g)
      else:
        l,o = self._recycle(p, backprop=False)
      _loss.append(l); _outs.append(o)
      _losses.append(o["losses"])

    # update
    if backprop:
      self._grad = jax.tree_map(lambda *v: jnp.stack(v).mean(0), *_grad)
    else:
      self._grad = jax.tree_map(jnp.zeros_like, self.params)

    self._loss = jnp.asarray(_loss).mean()
    self._losses = jax.tree_map(lambda *v: jnp.asarray(v).mean(), *_losses)
    self._outs = _outs[0]
    
    if callback is not None: callback(self)

  #-------------------------------------
  # STEP FUNCTION
  #-------------------------------------
  def _step(self, weights=None, opt=None,
            lr_scale=1.0, callback=None, **kwargs):
    '''do one step of gradient descent'''
    
    # update
    self.run(weights=weights, opt=opt, callback=callback)

    # normalize gradient
    g = self._grad["seq"]
    gn = jnp.linalg.norm(g,axis=(-1,-2),keepdims=True)
    self._grad["seq"] *= lr_scale * jnp.sqrt(self._len)/(gn+1e-7)

    # apply gradient
    self._state = self._update_fun(self._k, self._grad, self._state)
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

    if opt is not None: self.opt.update(opt)
    if weights is not None: self.opt["weights"].update(weights)
      
    self.set_opt(hard=hard, dropout=dropout, gumbel=gumbel)
    if e_soft is None: e_soft = soft
    if e_temp is None: e_temp = temp
    for i in range(iters):
      self.opt["temp"] = e_temp + (temp - e_temp) * (1-i/(iters-1)) ** 2
      self.opt["soft"] = soft + (e_soft - soft) * i/(iters-1)
      # decay learning rate based on temperature
      lr_scale = (1 - self.opt["soft"]) + (self.opt["soft"] * self.opt["temp"])
      self._step(lr_scale=lr_scale, callback=callback, **kwargs)

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
      self._step(**kwargs)
  
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

    for i in range(iters):
      seq = self._outs["seq"]["hard"]
      plddt = None
      if use_plddt:
        plddt = np.asarray(self._outs["plddt"])
        if self.protocol == "binder":
          plddt = plddt[self._target_len:]
        else:
          plddt = plddt[:self._len]
      
      buff = []
      for _ in range(tries):
        self.params["seq"] = mut(seq, plddt)
        self.run(backprop=False)
        buff.append({"outs":self._outs,"loss":self._loss,"losses":self._losses})
      
      # accept best
      buff = buff[np.argmin([x["loss"] for x in buff])]
      self._outs, self._loss, self._losses = buff["outs"], buff["loss"], buff["losses"]
      self._k += 1
      self._save_results(save_best=save_best, verbose=verbose)
