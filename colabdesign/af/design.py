import random, os
import jax
import jax.numpy as jnp
import numpy as np
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.shared.utils import update_dict, Key, dict_to_str

try:
  from jax.example_libraries.optimizers import sgd, adam
except:
  from jax.experimental.optimizers import sgd, adam

####################################################
# AF_DESIGN - design functions
####################################################
#\
# \_af_design
# |\
# | \_restart
# | |\_init_seq
# |  \_setup_optimizer
# |
#  \
#   \_design
#    \_step
#    |\_run
#    | \_recycle
#     \_save_results
#
####################################################

class _af_design:

  def restart(self, seed=None, weights=None, opt=None, set_defaults=False, keep_history=False,
              seq_init=None, mode=None, seq=None, rm_aa=None, add_seq=False,
              optimizer="sgd", **kwargs):    

    # set weights and options
    if set_defaults:
      update_dict(self._opt, opt)
      update_dict(self._opt["weights"], weights)
            
    if not keep_history:
      self.opt = jax.tree_map(lambda x:x, self._opt) # copy

    self.set_opt(opt)
    self.set_weights(weights)

    # initialize sequence
    self.key = Key(seed=seed).get
    self._init_seq(mode, seq, rm_aa, add_seq, seq_init)

    # setup optimizer
    self._setup_optimizer(optimizer, **kwargs)

    # initialize trajectory
    if not keep_history:
      self._traj = {"losses":[],"seq":[],"xyz":[],"plddt":[],"pae":[]}
      self.clear_best()
    
    # clear previous results
    if hasattr(self,"aux"): del self.aux  
    if hasattr(self,"loss"): del self.loss

  def _init_seq(self, mode=None, seq=None, rm_aa=None, add_seq=False, seq_init=None):
    '''initialize sequence'''

    # backward compatibility
    if seq_init is not None:
      print("WARNING: 'seq_init' is being deprecated, use options: 'mode' and 'seq'")
      x = seq_init
      if isinstance(x, (np.ndarray, jnp.ndarray)): seq = x
      elif isinstance(x, str):
        if len(x) == self._len: seq = x
        else: mode = x

    if mode is None: mode = []
    
    # initialize sequence
    shape = (self._num, self._len, 20)
    x = 0.01 * jax.random.normal(self.key(), shape)

    # initialize bias
    bias = jnp.zeros(shape[1:])
    use_bias = False

    # initialize with input sequence
    if "wildtype" in mode or "wt" in mode:
      seq = jax.nn.one_hot(self._wt_aatype, 20)
      if self.protocol == "partial":
        seq = jnp.zeros(shape).at[...,self.opt["pos"],:].set(seq)
      else:
        seq = seq[...,:self._len]
    if seq is not None:
      if isinstance(seq, str):
        seq = jnp.array([residue_constants.restype_order.get(aa,-1) for aa in seq])
        seq = jax.nn.one_hot(seq,20)
      x = x + seq
      if add_seq:
        bias = bias + seq * 1e8
        use_bias = True 

    if "gumbel" in mode: x += jax.random.gumbel(self.key(), shape) / self.opt["alpha"]
    if "soft" in mode:   x = jax.nn.softmax(x * self.opt["alpha"])

    # disable certain amino acids
    if rm_aa is not None:
      for aa in rm_aa.split(","):
        bias = bias - 1e8 * np.eye(20)[residue_constants.restype_order[aa]]
        use_bias = True

    # set bias and params
    if use_bias: self.opt["bias"] = bias
    self.params = {"seq":x}

  def _setup_optimizer(self, optimizer="sgd", **kwargs):
    '''setup which optimizer to use'''
    
    if optimizer == "adam":
      optimizer = adam
      self._opt["lr"] = 0.02
    
    if optimizer == "sgd":
      optimizer = sgd
      self._opt["lr"] = 0.1
      
    self._init_fun, self._update_fun, self._get_params = optimizer(1.0, **kwargs)
    self._state = self._init_fun(self.params)
    self._k = 0
  
  def design(self, iters,
             temp=1, e_temp=None,
             soft=0, e_soft=None,
             hard=0, e_hard=None,
             dropout=True, opt=None, weights=None,
             callback=None, save_best=False, verbose=1):
      
    # override settings if defined
    self.set_opt(opt)
    self.set_weights(weights)
    self.set_opt(dropout=dropout, hard=hard)

    if e_soft is None: e_soft = soft
    if e_hard is None: e_hard = hard
    if e_temp is None: e_temp = temp
    for i in range(iters):
      self.set_opt(soft=(soft+(e_soft-soft)*((i+1)/iters)),
                   hard=(hard+(e_hard-hard)*((i+1)/iters)),
                   temp=(e_temp+(temp-e_temp)*(1-(i+1)/iters)**2))
      
      # decay learning rate based on temperature
      lr_scale = (1 - self.opt["soft"]) + (self.opt["soft"] * self.opt["temp"])
      
      self.step(lr_scale=lr_scale, callback=callback,
                save_best=save_best, verbose=verbose)

  def step(self, lr_scale=1.0, backprop=True, callback=None,
           save_best=False, verbose=1):
    '''do one step of gradient descent'''
    
    # run
    self.run(backprop=backprop)
    if callback is not None: callback(self)

    # normalize gradient
    g = self.grad["seq"]
    gn = jnp.linalg.norm(g,axis=(-1,-2),keepdims=True)
    self.grad["seq"] *= jnp.sqrt(self._len)/(gn+1e-7)

    # set learning rate
    lr = self.opt["lr"] * lr_scale
    self.grad = jax.tree_map(lambda x:x*lr, self.grad)

    # apply gradient
    self._state = self._update_fun(self._k, self.grad, self._state)
    self.params = self._get_params(self._state)

    # increment
    self._k += 1

    # save results
    self.save_results(save_best=save_best, verbose=verbose)
    
  def run(self, backprop=True):
    '''run model to get outputs, losses and gradients'''
    
    # decide which model params to use
    m = self.opt["models"]
    ns = jnp.arange(2) if self._args["use_templates"] else jnp.arange(5)
    if self._args["sample_models"] and m != len(ns):
      model_num = jax.random.choice(self.key(),ns,(m,),replace=False)
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
      self.grad = jax.tree_map(lambda *x: jnp.stack(x).mean(0), *_grad)
    else:
      self.grad = jax.tree_map(jnp.zeros_like, self.params)

    # update loss (take mean across models)
    self.loss = jnp.stack(_loss).mean()
    
    # update aux
    _aux = jax.tree_map(lambda *x: jnp.stack(x), *_aux)    
    self.aux = jax.tree_map(lambda x:x[0], _aux)
    self.aux["losses"] = jax.tree_map(lambda x: x.mean(0), _aux["losses"])
    self.aux["model_num"] = model_num    
    self._outs = self.aux # backward compatibility
      
  def _recycle(self, model_params, backprop=True):                
    # initialize previous (for recycle)
    L = self._inputs["residue_index"].shape[-1]
    self._inputs["prev"] = {'prev_msa_first_row': np.zeros([L,256]),
                            'prev_pair': np.zeros([L,L,128]),
                            'prev_pos': np.zeros([L,37,3])}
              
    if self._args["recycle_mode"] in ["backprop","add_prev"]:
      self.opt["recycles"] = self._runner.config.model.num_recycle

    recycles = self.opt["recycles"]

    if self._args["recycle_mode"] == "average":
      # average gradient across recycles
      grad = []
      for r in range(recycles+1):
        if backprop:
          (loss, aux), _grad = self._grad_fn(self.params, model_params, self._inputs, self.key(), self.opt)
          grad.append(_grad)
        else:
          loss, aux = self._fn(self.params, model_params, self._inputs, self.key(), self.opt)
        self._inputs["prev"] = aux["prev"]
      if backprop:
        grad = jax.tree_map(lambda *x: jnp.asarray(x).mean(0), *grad)

    else:
      if self._args["recycle_mode"] == "sample":
        self.opt["recycles"] = jax.random.randint(self.key(), [], 0, recycles+1)
      
      if backprop:
        (loss, aux), grad = self._grad_fn(self.params, model_params, self._inputs, self.key(), self.opt)
      else:
        loss, aux = self._fn(self.params, model_params, self._inputs, self.key(), self.opt)

    aux["recycles"] = self.opt["recycles"]
    self.opt["recycles"] = recycles
    
    if backprop:
      return (loss, aux), grad
    else:
      return loss, aux
    
  def save_results(self, save_best=False, verbose=1):
    '''save the results and update trajectory'''

    # save best result
    if save_best and self.loss < self._best_loss:
      self._best_loss = self.loss
      self._best_aux = self._best_outs = self.aux

    # save losses
    losses = jax.tree_map(float,self.aux["losses"])

    losses.update({"models":         self.aux["model_num"],
                   "recycles":   int(self.aux["recycles"]),
                   "hard":     float(self.opt["hard"]),
                   "soft":     float(self.opt["soft"]),
                   "temp":     float(self.opt["temp"]),
                   "loss":     float(self.loss)})
    
    if self.protocol == "fixbb" or (self.protocol == "binder" and self._redesign):
      # compute sequence recovery
      _aatype = self.aux["seq"]["hard"].argmax(-1)
      L = min(_aatype.shape[-1], self._wt_aatype.shape[-1])
      losses["seqid"] = float((_aatype[...,:L] == self._wt_aatype[...,:L]).mean())

    # print results
    if verbose and (self._k % verbose) == 0:
      keys = ["models","recycles","hard","soft","temp","seqid","loss",
              "msa_ent","plddt","pae","helix","con","i_pae","i_con",
              "sc_fape","sc_rmsd","dgram_cce","fape","rmsd"]
      print(dict_to_str(losses, filt=self.opt["weights"], print_str=f"{self._k}",
                        keys=keys, ok="rmsd"))

    # save trajectory
    traj = {"losses":losses, "seq":np.asarray(self.aux["seq"]["pseudo"])}
    traj["xyz"] = np.asarray(self.aux["final_atom_positions"][:,1,:])
    traj["plddt"] = np.asarray(self.aux["plddt"])
    if "pae" in self.aux: traj["pae"] = np.asarray(self.aux["pae"])
    for k,v in traj.items(): self._traj[k].append(v)

  def clear_best(self):
    self._best_loss = np.inf
    self._best_aux = self._best_outs = None

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
                    temp=1.0, dropout=True, **kwargs):
    '''two stage design (soft→hard)'''
    self.design(soft_iters, soft=True, temp=temp, dropout=dropout, **kwargs)
    self.design(temp_iters, soft=True, temp=temp, dropout=dropout,  e_temp=1e-2, **kwargs)
    self.design(hard_iters, soft=True, temp=1e-2, dropout=False, hard=True, save_best=True, **kwargs)

  def design_3stage(self, soft_iters=300, temp_iters=100, hard_iters=10, waves=1,
                    temp=1.0, dropout=True, **kwargs):
    '''three stage design (logits→soft→hard)*waves'''
    M = 2 if self._args["use_templates"] else 5
    for w in range(waves):
      if hasattr(self,"aux"):
        self.params["seq"] = self.aux["seq"]["pseudo"]
        self._state = self._init_fun(self.params)
        self._k = 0
      self.set_opt(models=1) # sample models
      self.design(soft_iters, e_soft=True, temp=temp, dropout=dropout, **kwargs)
      self.design(temp_iters, soft=True,   temp=temp, dropout=dropout, e_temp=1e-2,**kwargs)
      self.set_opt(models=M) # use all models
      self.design(hard_iters, soft=True,   temp=1e-2, dropout=False, hard=True, save_best=True, **kwargs)

  ######################################################################################################
  # EXPERIMENTAL STUFF BELOW
  ######################################################################################################

  def design_semigreedy(self, iters=100, tries=20, dropout=False,
                        use_plddt=True, save_best=True, verbose=1):
    '''semigreedy search'''
    
    self.set_opt(hard=True, soft=True, temp=1.0, dropout=dropout)
    if self._k == 0:
      self.run(backprop=False)

    def mut(seq, plddt=None):
      '''mutate random position'''
      L,A = seq.shape[-2:]
      while True:
        if plddt is None:
          i = jax.random.randint(self.key(),[],0,L)
        else:
          p = (1-plddt)/(1-plddt).sum(-1,keepdims=True)
          # bias mutations towards positions with low pLDDT
          # https://www.biorxiv.org/content/10.1101/2021.08.24.457549v1
          i = jax.random.choice(self.key(),jnp.arange(L),[],p=p)

        a = jax.random.randint(self.key(),[],0,A)
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
      self.save_results(save_best=save_best, verbose=verbose)

  def template_predesign(self, iters=100, soft=True, hard=False, 
                         temp=1.0, dropout=True, **kwargs):
    '''use template for predesign stage, then increase template dropout until gone'''
    self.set_opt(hard=hard, soft=soft, dropout=dropout, temp=temp)
    for i in range(iters):
      self.opt["template"]["dropout"] = (i+1)/iters
      self.step(**kwargs)