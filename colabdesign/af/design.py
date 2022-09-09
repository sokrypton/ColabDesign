import random, os
import jax
import jax.numpy as jnp
import numpy as np
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.shared.utils import copy_dict, update_dict, Key, dict_to_str, to_float

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
#  \
#   \_design
#    \_step
#     \_run
#      \_recycle
#       \_single
#
####################################################

class _af_design:

  def restart(self, seed=None, optimizer="sgd", opt=None, weights=None,
              seq=None, keep_history=False, reset_opt=True, **kwargs):   
    '''
    restart the optimization
    ------------
    note: model.restart() resets the [opt]ions and weights to their defaults
    use model.set_opt(..., set_defaults=True) and model.set_weights(..., set_defaults=True)
    or model.restart(reset_opt=False) to avoid this
    ------------
    seed=0 - set seed for reproducibility
    reset_opt=False - do NOT reset [opt]ions/weights to defaults
    keep_history=True - do NOT clear the trajectory/[opt]ions/weights
    '''
    # reset [opt]ions
    if reset_opt and not keep_history:
      self.opt = copy_dict(self._opt)

    # update options/settings (if defined)
    self.set_opt(opt)
    self.set_weights(weights)

    # setup optimizer
    if isinstance(optimizer, str):
      if optimizer == "adam": (optimizer,lr) = (adam,0.02)
      if optimizer == "sgd":  (optimizer,lr) = (sgd,0.1)
      self.opt["lr"] = lr

    self._init_fun, self._update_fun, self._get_params = optimizer(1.0)

    # initialize sequence
    self.key = Key(seed=seed).get
    self.set_seq(seq=seq, **kwargs)
    self._k = 0

    if not keep_history:
      # initialize trajectory
      self._traj = {"log":[],"seq":[],"xyz":[],"plddt":[],"pae":[]}
      self._best, self._tmp = {}, {}
      
  def run(self, num_recycles=None, backprop=True, callback=None):
    '''run model to get outputs, losses and gradients'''

    callbacks = [callback]
    if self._args["use_crop"]: callbacks.append(self._crop())
    
    # decide which model params to use
    ns_name = self._model_names.copy()
    ns = list(range(len(ns_name)))

    # sub select number of model params
    if self._args["models"] is not None:
      models = self._args["models"]
      models = [models] if isinstance(models,(int,str)) else models
      ns = [ns[n if isinstance(n,int) else ns_name.index(n)] for n in models]
    
    ns = jnp.array(ns)
    m = min(self.opt["num_models"],len(ns))
    if self.opt["sample_models"] and m != len(ns):
      model_num = jax.random.choice(self.key(),ns,(m,),replace=False)
    else:
      model_num = ns[:m]      
    model_num = np.array(model_num).tolist()

    # loop through model params
    aux = []
    for n in model_num:
      p = self._model_params[n]
      aux.append(self._recycle(p, num_recycles=num_recycles, backprop=backprop))
    aux = jax.tree_map(lambda *x: jnp.stack(x), *aux)

    # update aux
    self.aux = jax.tree_map(lambda x:x[0], aux)
    self.aux["all"] = aux

    # average losses and gradients
    self.aux["loss"] = aux["loss"].mean()
    self.aux["losses"] = jax.tree_map(lambda x: x.mean(0), aux["losses"])
    self.aux["grad"] = jax.tree_map(lambda x: x.mean(0), aux["grad"])

    # callback
    for callback in callbacks:
      if callback is not None: callback(self)

    # update log
    self.aux["log"] = {**self.aux["losses"], "loss":self.aux["loss"],
                       "ptm":self.aux["ptm"], "i_ptm":self.aux["i_ptm"]}
    self.aux["log"].update({k:self.opt[k] for k in ["hard","soft","temp"]})

    # compute sequence recovery
    if self.protocol in ["fixbb","partial"] or (self.protocol == "binder" and self._args["redesign"]):
      if self.protocol == "partial":
        aatype = self.aux["aatype"].argmax(-1)[...,self.opt["pos"]]
      else:
        aatype = self.aux["seq"]["pseudo"].argmax(-1)

      mask = self._wt_aatype != -1
      true = self._wt_aatype[mask]
      pred = aatype[...,mask]
      self.aux["log"]["seqid"] = (true == pred).mean()

    self.aux["log"] = to_float(self.aux["log"])
    self.aux["log"].update({"recycles":int(self.aux["num_recycles"]),
                            "models":model_num})

  def _single(self, model_params, backprop=True):
    '''single pass through the model'''
    flags  = [self._params, model_params, self._inputs, self.key(), self.opt]
    if backprop:
      (loss, aux), grad = self._model["grad_fn"](*flags)
    else:
      loss, aux = self._model["fn"](*flags)
      grad = jax.tree_map(jnp.zeros_like, self._params)
    aux.update({"loss":loss,"grad":grad})
    return aux

  def _recycle(self, model_params, num_recycles=None, backprop=True):   
    '''multiple passes through the model (aka recycle)'''
    mode = self._args["recycle_mode"]
    if num_recycles is None:
      num_recycles = self.opt["num_recycles"]

    if mode in ["backprop","add_prev"]:
      # recycles compiled into model, only need single-pass
      aux = self._single(model_params, backprop)
    
    else:
      L = self._inputs["residue_index"].shape[0]
      if self._args["use_crop"]: L = self.opt["crop_pos"].shape[0]
      
      # intialize previous
      self._inputs["prev"] = {'prev_msa_first_row': np.zeros([L,256]),
                              'prev_pair': np.zeros([L,L,128]),
                              'prev_pos': np.zeros([L,37,3])}

      # decide which layers to compute gradients for
      cycles = (num_recycles + 1)
      mask = [0] * cycles
      if mode == "sample":  mask[jax.random.randint(self.key(),[],0,cycles)] = 1
      if mode == "average": mask = [1/cycles] * cycles
      if mode == "last":    mask[-1] = 1
      if mode == "first":   mask[0] = 1
      
      # gather gradients across recycles 
      grad = []
      for m in mask:
        if m == 0:
          aux = self._single(model_params, backprop=False)
        else:
          aux = self._single(model_params, backprop)
          grad.append(jax.tree_map(lambda x:x*m, aux["grad"]))
        self._inputs["prev"] = aux["prev"]

      aux["grad"] = jax.tree_map(lambda *x: jnp.stack(x).sum(0), *grad)
    
    aux["num_recycles"] = num_recycles
    return aux

  def step(self, lr_scale=1.0, num_recycles=None, backprop=True,
           callback=None, stats_correct=False, save_best=False, verbose=1):
    '''do one step of gradient descent'''
    
    # run
    self.run(num_recycles=num_recycles, backprop=backprop, callback=callback)

    # apply gradient
    g = self.aux["grad"]["seq"]
    eff_len = (jnp.square(g).sum(-1,keepdims=True) > 0).sum(-2,keepdims=True)
    
    # statistical correction - doi:10.1101/2022.04.29.490102
    if stats_correct:
      g = g - g.sum(-2,keepdims=True) / eff_len

    # normalize gradient
    gn = jnp.linalg.norm(g,axis=(-1,-2),keepdims=True)
    self.aux["grad"]["seq"] = g * jnp.sqrt(eff_len)/(gn+1e-7)

    # set learning rate
    lr = self.opt["lr"] * lr_scale
    self.aux["grad"] = jax.tree_map(lambda x:x*lr, self.aux["grad"])

    # update state/params
    self._state = self._update_fun(self._k, self.aux["grad"], self._state)
    self._params = self._get_params(self._state)

    # increment
    self._k += 1

    # save results
    self._save_results(save_best=save_best, verbose=verbose)

  def _update_traj(self):
    traj = {"seq":   self.aux["seq"]["pseudo"],
            "xyz":   self.aux["atom_positions"][:,1,:],
            "plddt": self.aux["plddt"],
            "pae":   self.aux["pae"]}
    traj = jax.tree_map(np.array, traj)
    traj["log"] = self.aux["log"]
    for k,v in traj.items(): self._traj[k].append(v)

  def _print_log(self, print_str=None):
    keys = ["models","recycles","hard","soft","temp","seqid","loss",
            "seq_ent","mlm","plddt","pae","exp_res","con","i_con",
            "sc_fape","sc_rmsd","dgram_cce","fape","ptm"]
    if sum(self._lengths) > 1: keys.append("i_ptm")
    keys.append("rmsd")
    print(dict_to_str(self.aux["log"], filt=self.opt["weights"],
                      print_str=print_str, keys=keys, ok="rmsd"))

  def _save_best(self):
    metric = self.aux["log"][self._args["best_metric"]]
    if "metric" not in self._best or metric < self._best["metric"]:
      self._best.update({"metric":metric, "aux":self.aux})

  def _save_results(self, save_best=False, verbose=True):    
    self._update_traj()
    if save_best: self._save_best()
    if verbose and (self._k % verbose) == 0:
      self._print_log(f"{self._k}")

  def predict(self, seq=None, num_recycles=None, num_models=None, models=None, verbose=True):
    '''predict structure for input sequence (if provided)'''
    
    # save settings
    (opt, args, params) = (copy_dict(x) for x in [self.opt, self._args, self._params])    

    # set settings
    if seq is not None: self.set_seq(seq=seq, set_state=False)
    if models is None:
      models = self._model_names if num_models is None else self._model_names[:num_models]
    num_models = len(models) if isinstance(models,list) else 1
    
    self.set_opt(hard=True, dropout=False, sample_models=False,
                 models=models, num_models=num_models,
                 mlm_dropout=0.0, use_crop=False)
    # run
    self.run(num_recycles=num_recycles, backprop=False)
    if verbose: self._print_log("predict")

    # reset settings
    (self.opt, self._args, self._params) = (opt, args, params)

  # ---------------------------------------------------------------------------------
  # example design functions
  # ---------------------------------------------------------------------------------
  def design(self, iters=100,
             soft=0.0, e_soft=None,
             temp=1.0, e_temp=None,
             hard=0.0, e_hard=None,
             step=1.0, e_step=None,
             dropout=True, opt=None, weights=None,
             mlm_dropout=0.05, num_recycles=None,
             backprop=True, callback=None,
             save_best=False, verbose=1):

    # update options/settings (if defined)
    self.set_opt(opt, dropout=dropout, mlm_dropout=mlm_dropout)
    self.set_weights(weights)
    
    m = {"soft":[soft,e_soft],"temp":[temp,e_temp],
         "hard":[hard,e_hard],"step":[step,e_step]}
    m = {k:[s,(s if e is None else e)] for k,(s,e) in m.items()}

    for i in range(iters):
      for k,(s,e) in m.items():
        if k == "temp":
          self.set_opt({k:(e+(s-e)*(1-(i+1)/iters)**2)})
        else:
          v = (s+(e-s)*((i+1)/iters))
          if k == "step": step = v
          else: self.set_opt({k:v})
      
      # decay learning rate based on temperature
      lr_scale = step * ((1 - self.opt["soft"]) + (self.opt["soft"] * self.opt["temp"]))
      
      self.step(lr_scale=lr_scale, num_recycles=num_recycles, backprop=backprop,
                callback=callback, save_best=save_best, verbose=verbose)

  def design_logits(self, iters=100, **kwargs):
    ''' optimize logits '''
    self.design(iters, **kwargs)

  def design_soft(self, iters=100, temp=1, **kwargs):
    ''' optimize softmax(logits/temp)'''
    self.design(iters, soft=1, temp=temp, **kwargs)
  
  def design_hard(self, iters=100, **kwargs):
    ''' optimize argmax(logits) '''
    self.design(iters, soft=1, hard=1, **kwargs)

  # ---------------------------------------------------------------------------------
  # experimental
  # ---------------------------------------------------------------------------------
  def design_3stage(self, soft_iters=300, temp_iters=100, hard_iters=10,
                    ramp_recycles=True, num_recycles=None, num_models=1, **kwargs):
    '''three stage design (logits→soft→hard)'''

    # set starting options
    kwargs["num_recycles"] = num_recycles
    self.set_opt(num_models=num_models, sample_models=True)

    # logits -> softmax(logits/1.0)
    if ramp_recycles and self._args["recycle_mode"] not in ["add_prev","backprop"]:
      R = self.opt["num_recycles"] if num_recycles is None else num_recycles
      p = 1.0 / (R + 1)
      iters = soft_iters // (R + 1)
      for r in range(R + 1):
        kwargs["num_recycles"] = r
        self.design_logits(iters, soft=r*p, e_soft=(r+1)*p, **kwargs)    
    else:
      self.design_logits(soft_iters, e_soft=1, **kwargs)
    
    # softmax(logits/1.0) -> softmax(logits/0.01)
    self.design_soft(temp_iters, e_temp=1e-2, **kwargs)

    self.set_opt(num_models=len(self._model_params)) # use all models
    self.design_hard(hard_iters, temp=1e-2, dropout=False, mlm_dropout=0.0, save_best=True,
                     num_recycles=kwargs["num_recycles"], verbose=kwargs.get("verbose",1))

  def design_semigreedy(self, iters=100, tries=20, num_recycles=None, num_models=1,
                        use_plddt=True, save_best=True, verbose=1):
    
    '''semigreedy search'''    
    self.set_opt(hard=True, dropout=False, use_crop=False,
                 num_models=num_models, sample_models=False)
    
    if self._k == 0:
      self.run(num_recycles=num_recycles, backprop=False)

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
      return jax.nn.one_hot(self._params["seq"].argmax(-1),20)

    # optimize!
    for i in range(iters):
      seq = get_seq()
      plddt = None
      if use_plddt:
        plddt = np.asarray(self.aux["all"]["plddt"].mean(0))
        if self.protocol == "binder":
          plddt = plddt[self._target_len:]
        else:
          plddt = plddt[:self._len]
      
      buff = []
      for _ in range(tries):
        self.set_seq(seq=mut(seq, plddt), set_state=False)
        self.run(num_recycles=num_recycles, backprop=False)
        buff.append({"aux":self.aux, "seq":self._params["seq"]})
      
      # accept best      
      buff = buff[np.argmin([x["aux"]["loss"] for x in buff])]
      
      self.set_seq(seq=buff["seq"])
      self.aux = buff["aux"]
      self._k += 1
      self._save_results(save_best=save_best, verbose=verbose)

  # ---------------------------------------------------------------------------------

  def template_predesign(self, iters=100, soft=True, hard=False, 
                         temp=1.0, dropout=True, **kwargs):
    '''use template for predesign stage, then increase template dropout until gone'''
    self.set_opt(hard=hard, soft=soft, dropout=dropout, temp=temp)
    for i in range(iters):
      self.opt["template"]["dropout"] = (i+1)/iters
      self.step(**kwargs)
