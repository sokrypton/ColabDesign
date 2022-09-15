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
    self.set_seed(seed)
    self.set_seq(seq=seq, **kwargs)
    self._k = 0

    if not keep_history:
      # initialize trajectory
      self._traj = {"log":[],"seq":[],"xyz":[],"plddt":[],"pae":[]}
      self._best, self._tmp = {}, {}

  def _get_model_nums(self, num_models=None, sample_models=None, models=None):
    '''decide which model params to use'''
    if num_models is None: num_models = self.opt["num_models"]
    if sample_models is None: sample_models = self.opt["sample_models"]

    ns_name = self._model_names
    ns = list(range(len(ns_name)))
    if models is not None:
      models = models if isinstance(models,list) else [models]
      ns = [ns[n if isinstance(n,int) else ns_name.index(n)] for n in models]
    
    m = min(num_models,len(ns))
    if sample_models and m != len(ns):
      ns = jnp.array(ns)
      model_nums = jax.random.choice(self.key(),ns,(m,),replace=False)
      model_nums = np.array(model_nums).tolist()
    else:
      model_nums = ns[:m]

    return model_nums    

  def run(self, num_recycles=None, num_models=None, sample_models=None, models=None,
          backprop=True, callback=None, model_nums=None, return_aux=False):
    '''run model to get outputs, losses and gradients'''

    callbacks = [callback]
    if self._args["use_crop"]: callbacks.append(self._crop())
    
    # decide which model params to use
    if model_nums is None:
      model_nums = self._get_model_nums(num_models, sample_models, models)
  
    # loop through model params
    auxs = []
    for n in model_nums:
      p = self._model_params[n]
      auxs.append(self._recycle(p, num_recycles=num_recycles, backprop=backprop))
    auxs = jax.tree_map(lambda *x: np.stack(x), *auxs)

    # update aux
    aux = jax.tree_map(lambda x:x[0], auxs)
    aux["all"] = auxs

    # average losses and gradients
    aux["loss"] = auxs["loss"].mean()
    aux["losses"] = jax.tree_map(lambda x: x.mean(0), auxs["losses"])
    aux["grad"] = jax.tree_map(lambda x: x.mean(0), auxs["grad"])

    # callback
    for callback in callbacks:
      if callback is not None: callback(self)

    # update log
    aux["log"] = {**aux["losses"], "loss":aux["loss"]}
    aux["log"].update({k:auxs[k].mean() for k in ["ptm","i_ptm"]})
    aux["log"].update({k:self.opt[k] for k in ["hard","soft","temp"]})
    aux["log"]["plddt"] = 1 - aux["log"]["plddt"]

    # compute sequence recovery
    if self.protocol in ["fixbb","partial"] or (self.protocol == "binder" and self._args["redesign"]):
      if self.protocol == "partial":
        aatype = aux["aatype"].argmax(-1)[...,self.opt["pos"]]
      else:
        aatype = aux["seq"]["pseudo"].argmax(-1)

      mask = self._wt_aatype != -1
      true = self._wt_aatype[mask]
      pred = aatype[...,mask]
      aux["log"]["seqid"] = (true == pred).mean()

    aux["log"] = to_float(aux["log"])
    aux["log"].update({"recycles":int(aux["num_recycles"]), "models":model_nums})
    
    if return_aux:
      return aux
    else:
      self.aux = aux

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

  def step(self, lr_scale=1.0, num_recycles=None,
           num_models=None, sample_models=None, models=None, backprop=True,
           callback=None, stats_correct=False, save_best=False, verbose=1):
    '''do one step of gradient descent'''
    
    # run
    self.run(num_recycles=num_recycles, num_models=num_models, sample_models=sample_models,
             models=models, backprop=backprop, callback=callback)

    # apply gradient
    g = self.aux["grad"]["seq"]
    eff_len = (np.square(g).sum(-1,keepdims=True) > 0).sum(-2,keepdims=True)
    
    # statistical correction - doi:10.1101/2022.04.29.490102
    if stats_correct:
      g = g - g.sum(-2,keepdims=True) / eff_len

    # normalize gradient
    gn = np.linalg.norm(g,axis=(-1,-2),keepdims=True)
    self.aux["grad"]["seq"] = g * np.sqrt(eff_len)/(gn+1e-7)

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

  def _print_log(self, print_str=None, aux=None):
    if aux is None: aux = self.aux
    keys = ["models","recycles","hard","soft","temp","seqid","loss",
            "seq_ent","mlm","pae","exp_res","con","i_con",
            "sc_fape","sc_rmsd","dgram_cce","fape","plddt","ptm"]
    
    if "i_ptm" in aux["log"]:
      if len(self._lengths) > 1:
        keys.append("i_ptm")
      else:
        aux["log"].pop("i_ptm")

    print(dict_to_str(aux["log"], filt=self.opt["weights"],
                      print_str=print_str, keys=keys+["rmsd"], ok=["plddt","rmsd"]))

  def _save_results(self, aux=None, save_best=False, verbose=True):
    if aux is None: aux = self.aux
    # update traj
    traj = {"seq":   aux["seq"]["pseudo"],
            "xyz":   aux["atom_positions"][:,1,:],
            "plddt": aux["plddt"],
            "pae":   aux["pae"]}
    traj = jax.tree_map(np.array, traj)
    traj["log"] = aux["log"]
    for k,v in traj.items(): self._traj[k].append(v)

    # save best
    if save_best:
      metric = aux["log"][self._args["best_metric"]]
      if "metric" not in self._best or metric < self._best["metric"]:
        self._best["aux"] = aux
        self._best["metric"] = metric

    if verbose:
      self._print_log(f"{self._k}", aux=aux)

  def predict(self, seq=None,
              num_models=None, num_recycles=None, models=None, sample_models=False,
              dropout=False, mlm_dropout=False, hard=True, soft=False, temp=1,
              return_aux=False, verbose=True,  seed=None, **kwargs):
    '''predict structure for input sequence (if provided)'''

    # set seed if defined
    if seed is not None: self.set_seed(seed)

    # save settings
    (opt, args, params) = (copy_dict(x) for x in [self.opt, self._args, self._params])

    # set [seq]uence/[opt]ions
    if seq is not None: self.set_seq(seq=seq, set_state=False)    
    self.set_opt(hard=hard, soft=soft, temp=temp, dropout=dropout,
                 mlm_dropout=mlm_dropout, use_crop=False, use_pssm=False)
        
    # run
    aux = self.run(num_recycles=num_recycles, num_models=num_models,
                   sample_models=sample_models, models=models, backprop=False,
                   return_aux=True, **kwargs)
    if verbose: self._print_log("predict", aux=aux)

    # reset settings
    (self.opt, self._args, self._params) = (opt, args, params)

    # return (or save) results
    if return_aux: return aux
    else: self.aux = aux

  # ---------------------------------------------------------------------------------
  # example design functions
  # ---------------------------------------------------------------------------------
  def design(self, iters=100,
             soft=0.0, e_soft=None,
             temp=1.0, e_temp=None,
             hard=0.0, e_hard=None,
             step=1.0, e_step=None,
             dropout=True, opt=None, weights=None, mlm_dropout=0.05, 
             num_recycles=None, num_models=None, sample_models=None, models=None,
             backprop=True, callback=None, save_best=False, verbose=1):

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
      
      self.step(lr_scale=lr_scale, num_recycles=num_recycles,
                num_models=num_models, sample_models=sample_models, models=models,
                backprop=backprop, callback=callback, save_best=save_best, verbose=verbose)

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
                    ramp_recycles=True, **kwargs):
    '''three stage design (logits→soft→hard)'''

    verbose = kwargs.get("verbose",1)
    if soft_iters > 0:
      if verbose: print("stage1: running logits→soft")
      # stage 1: logits -> softmax(logits/1.0)
      if ramp_recycles and self._args["recycle_mode"] not in ["add_prev","backprop"]:
        num_cycles = kwargs.pop("num_recycles", self.opt["num_recycles"]) + 1
        p = 1.0 / num_cycles
        for c in range(num_cycles):
          if verbose > 1 and c > 0: print(f"increasing number of recycles to {c+1}")
          kwargs["num_recycles"] = c
          iters = soft_iters // num_cycles
          self.design_logits(iters, soft=c*p, e_soft=(c+1)*p, **kwargs)        
      else:
        self.design_logits(soft_iters, e_soft=1, **kwargs)
    
    self._tmp["seq_logits"] = self.aux["seq"]["logits"]
    
    if temp_iters > 0:
      if verbose: print("stage2: running softmax(logits/1.0)→soft(logits/0.01)")
      # stage 2: softmax(logits/1.0) -> softmax(logits/0.01)
      self.design_soft(temp_iters, e_temp=1e-2, **kwargs)
    
    if hard_iters > 0:
      # stage 3: 
      if verbose: print("stage 3: running one_hot(logits)")
      for k in ["dropout","mlm_dropout","save_best","num_models"]: kwargs.pop(k,None)
      self.design_hard(hard_iters, temp=1e-2, dropout=False, mlm_dropout=0.0, save_best=True,
                      num_models=len(self._model_names), **kwargs)

  def design_semigreedy(self, iters=100, tries=20, dropout=False, use_plddt=True,
                        save_best=True, seq_logits=None, **kwargs):
    '''semigreedy search'''    

    def mut(seq, plddt=None, logits=None):
      '''mutate random position'''
      N,L,A = seq.shape

      # fix some positions
      i_prob = jnp.ones(L) if plddt is None else jax.nn.relu(1-plddt)
      if "fix_pos" in self.opt:
        seq, p = self._fix_pos(seq, return_p=True)
        i_prob = i_prob.at[p].set(0)
      
      # sample position
      # https://www.biorxiv.org/content/10.1101/2021.08.24.457549v1
      i = jax.random.choice(self.key(),jnp.arange(L),[],p=i_prob/i_prob.sum())

      # sample amino acid
      logits = np.array(0 if logits is None else logits)
      if logits.ndim == 3:   b = logits[:,i]
      elif logits.ndim == 2: b = logits[i]
      else:                  b = logits
      a_logits = b - seq[:,i] * 1e8
      a = jax.random.categorical(self.key(),a_logits)

      # return mutant
      return seq.at[:,i,:].set(jax.nn.one_hot(a,A))

    verbose = kwargs.pop("verbose",1)
    plddt = None
    seq = jax.nn.one_hot(self._params["seq"].argmax(-1),20)
    model_flags = {k:kwargs.pop(k,None) for k in ["num_models","sample_models","models"]}

    # optimize!
    for i in range(iters):

      buff = []
      model_nums = self._get_model_nums(**model_flags)
      for t in range(tries):
        mut_seq = mut(seq, plddt, logits=seq_logits)
        aux = self.predict(mut_seq, return_aux=True, model_nums=model_nums,
                           verbose=verbose>1, **kwargs)
        buff.append({"aux":aux, "seq":np.array(mut_seq)})

      # accept best
      losses = [x["aux"]["loss"] for x in buff]
      best = buff[np.argmin(losses)]
      self.aux, seq, self._k = best["aux"], jnp.array(best["seq"]), self._k + 1
      self.set_seq(seq=seq)
      self._save_results(save_best=save_best, verbose=verbose)

      # update plddt
      if use_plddt:
        plddt = best["aux"]["all"]["plddt"].mean(0)
        plddt = plddt[self._target_len:] if self.protocol == "binder" else plddt[:self._len]

  def design_pssm_semigreedy(self, soft_iters=300, hard_iters=32, tries=10,
                             ramp_recycles=True, ramp_models=True, **kwargs):

    verbose = kwargs.get("verbose",1)
    # stage 1: logits -> softmax(logits)
    if soft_iters > 0:
      self.design_3stage(soft_iters, 0, 0, ramp_recycles=ramp_recycles, **kwargs)
      self._tmp["seq_logits"] = kwargs["seq_logits"] = self.aux["seq"]["logits"]

    # stage 2: semi_greedy
    if hard_iters > 0:
      if verbose: print("running semigreedy optimization")
      kwargs["dropout"] = False
      if ramp_models:
        num_models = len(kwargs.get("models",self._model_names))
        iters = hard_iters
        for m in range(num_models):
          if verbose > 1 and m > 0: print(f'increasing number of models to {m+1}')
          kwargs["num_models"] = m + 1
          kwargs["save_best"] = (m + 1) == num_models
          self.design_semigreedy(iters, tries=tries, **kwargs)
          if m < 2: iters = iters // 2
      else:
        self.design_semigreedy(hard_iters, tries=tries, **kwargs)

  # ---------------------------------------------------------------------------------