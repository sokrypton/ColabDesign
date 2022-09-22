import random, os
import jax
import jax.numpy as jnp
import numpy as np
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.shared.utils import copy_dict, update_dict, Key, dict_to_str, to_float, softmax, categorical

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

  def restart(self, seed=None, opt=None, weights=None,
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
    if not keep_history:
      # initialize trajectory
      self._traj = {"log":[],"seq":[],"xyz":[],"plddt":[],"pae":[]}
      self._best, self._tmp = {}, {}

    # update options/settings (if defined)
    self.set_opt(opt)
    self.set_weights(weights)
  
    # initialize sequence
    self.set_seed(seed)
    self.set_seq(seq=seq, **kwargs)

    # reset optimizer
    self._k = 0
    self.set_optimizer()

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
      ns = np.array(ns)

      model_nums = np.random.choice(ns,(m,),replace=False)
      model_nums = np.array(model_nums).tolist()
    else:
      model_nums = ns[:m]

    return model_nums    

  def run(self, num_recycles=None, num_models=None, sample_models=None, models=None,
          backprop=True, callback=None, model_nums=None, return_aux=False):
    '''run model to get outputs, losses and gradients'''
    
    # decide which model params to use
    if model_nums is None:
      model_nums = self._get_model_nums(num_models, sample_models, models)
    assert len(model_nums) > 0, "ERROR: no model params defined"

    # loop through model params
    auxs = []
    for n in model_nums:
      p = self._model_params[n]
      auxs.append(self._recycle(p, num_recycles=num_recycles, backprop=backprop))
    auxs = jax.tree_map(lambda *x: np.stack(x), *auxs)

    # update aux
    self.aux = jax.tree_map(lambda x:x[0], auxs)
    self.aux["all"] = auxs

    # average losses and gradients
    self.aux["loss"] = auxs["loss"].mean()
    self.aux["losses"] = jax.tree_map(lambda x: x.mean(0), auxs["losses"])
    self.aux["grad"] = jax.tree_map(lambda x: x.mean(0), auxs["grad"])
    
    # design callbacks
    def append(x,fns):
      if not isinstance(fns,list): fns = [fns]
      x += [fn for fn in fns if fn is not None]
    callbacks = []
    append(callbacks, callback)
    append(callbacks, self._callbacks["design"])
    if self._args["use_crop"]: append(callbacks, self._crop())
    for callback in callbacks: callback(self)

    # update log
    self.aux["log"] = {**self.aux["losses"], "loss":self.aux["loss"]}
    self.aux["log"].update({k:auxs[k].mean() for k in ["ptm","i_ptm"]})
    self.aux["log"].update({k:self.opt[k] for k in ["hard","soft","temp"]})
    self.aux["log"]["plddt"] = 1 - self.aux["log"]["plddt"]

    # compute sequence recovery
    if self.protocol in ["fixbb","partial"] or (self.protocol == "binder" and self._args["redesign"]):
      if self.protocol == "partial":
        aatype = self.aux["aatype"][...,self.opt["pos"]]
      else:
        aatype = self.aux["seq"]["pseudo"].argmax(-1)

      mask = self._wt_aatype != -1
      true = self._wt_aatype[mask]
      pred = aatype[...,mask]
      self.aux["log"]["seqid"] = (true == pred).mean()

    self.aux["log"] = to_float(self.aux["log"])
    self.aux["log"].update({"recycles":int(self.aux["num_recycles"]), "models":model_nums})
    
    if return_aux: return self.aux

  def _single(self, model_params, backprop=True):
    '''single pass through the model'''
    self._inputs["opt"] = self.opt
    flags  = [self._params, model_params, self._inputs, self.key()]
    if backprop:
      (loss, aux), grad = self._model["grad_fn"](*flags)
    else:
      loss, aux = self._model["fn"](*flags)
      grad = jax.tree_map(np.zeros_like, self._params)
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

      if mode == "sample":  mask[np.random.randint(0,cycles)] = 1
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

      aux["grad"] = jax.tree_map(lambda *x: np.stack(x).sum(0), *grad)
    
    aux["num_recycles"] = num_recycles
    return aux

  def step(self, lr_scale=1.0, num_recycles=None,
           num_models=None, sample_models=None, models=None, backprop=True,
           callback=None, save_best=False, verbose=1):
    '''do one step of gradient descent'''
    
    # run
    self.run(num_recycles=num_recycles, num_models=num_models, sample_models=sample_models,
             models=models, backprop=backprop, callback=callback)

    # modify gradients    
    if self.opt["norm_seq_grad"]: self._norm_seq_grad()
    self._state, self.aux["grad"] = self._optimizer(self._state, self.aux["grad"], self._params)
  
    # apply gradients
    lr = self.opt["learning_rate"] * lr_scale
    self._params = jax.tree_map(lambda x,g:x-lr*g, self._params, self.aux["grad"])

    # increment
    self._k += 1

    # save results
    self._save_results(save_best=save_best, verbose=verbose)

  def _print_log(self, print_str=None, aux=None):
    if aux is None: aux = self.aux
    keys = ["models","recycles","hard","soft","temp","seqid","loss",
            "seq_ent","mlm","pae","i_pae","exp_res","con","i_con",
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

    if verbose and (self._k % verbose) == 0:
      self._print_log(f"{self._k}", aux=aux)

  def predict(self, seq=None,
              num_models=None, num_recycles=None, models=None, sample_models=False,
              dropout=False, hard=True, soft=False, temp=1,
              return_aux=False, verbose=True,  seed=None, **kwargs):
    '''predict structure for input sequence (if provided)'''

    # set seed if defined
    if seed is not None: self.set_seed(seed)

    # save settings
    (opt, args, params) = (copy_dict(x) for x in [self.opt, self._args, self._params])

    # set [seq]uence/[opt]ions
    if seq is not None: self.set_seq(seq=seq)    

    self.set_opt(hard=hard, soft=soft, temp=temp, dropout=dropout,
                 use_crop=False, use_pssm=False)
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
             dropout=True, opt=None, weights=None, 
             num_recycles=None, ramp_recycles=False, 
             num_models=None, sample_models=None, models=None,
             backprop=True, callback=None, save_best=False, verbose=1):

    # update options/settings (if defined)
    self.set_opt(opt, dropout=dropout)
    self.set_weights(weights)    
    m = {"soft":[soft,e_soft],"temp":[temp,e_temp],
         "hard":[hard,e_hard],"step":[step,e_step]}
    m = {k:[s,(s if e is None else e)] for k,(s,e) in m.items()}

    if ramp_recycles:
      if num_recycles is None:
        num_recycles = self.opt["num_recycles"]
      m["num_recycles"] = [0,num_recycles]

    for i in range(iters):
      for k,(s,e) in m.items():
        if k == "temp":
          self.set_opt({k:(e+(s-e)*(1-(i+1)/iters)**2)})
        else:
          v = (s+(e-s)*((i+1)/iters))
          if k == "step": step = v
          elif k == "num_recycles": num_recycles = round(v)
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

    # stage 1: logits -> softmax(logits/1.0)
    if soft_iters > 0:
      if verbose: print("Stage 1: running (logits → soft)")
      self.design_logits(soft_iters, e_soft=1,
        ramp_recycles=ramp_recycles, **kwargs)
      self._tmp["seq_logits"] = self.aux["seq"]["logits"]
      
    # stage 2: softmax(logits/1.0) -> softmax(logits/0.01)
    if temp_iters > 0:
      if verbose: print("Stage 2: running (soft → hard)")
      self.design_soft(temp_iters, e_temp=1e-2, **kwargs)
    
    # stage 3:
    if hard_iters > 0:
      if verbose: print("Stage 3: running (hard)")
      kwargs["dropout"] = False
      kwargs["save_best"] = True
      kwargs["num_models"] = len(self._model_names)
      self.design_hard(hard_iters, temp=1e-2, **kwargs)

  def _mutate(self, seq, plddt=None, logits=None, mutation_rate=1):
    '''mutate random position'''
    seq = np.array(seq)
    N,L = seq.shape

    # fix some positions
    i_prob = np.ones(L) if plddt is None else np.maximum(1-plddt,0)
    if "fix_pos" in self.opt:
      if "pos" in self.opt:
        p = self.opt["pos"][self.opt["fix_pos"]]
        seq[...,p] = self._wt_aatype_sub
      else:
        p = self.opt["fix_pos"]
        seq[...,p] = self._wt_aatype[...,p]
      i_prob[p] = 0
    
    for m in range(mutation_rate):
      # sample position
      # https://www.biorxiv.org/content/10.1101/2021.08.24.457549v1
      i = np.random.choice(np.arange(L),p=i_prob/i_prob.sum())

      # sample amino acid
      logits = np.array(0 if logits is None else logits)
      if logits.ndim == 3: logits = logits[:,i]
      elif logits.ndim == 2: logits = logits[i]
      a_logits = logits - np.eye(20)[seq[:,i]] * 1e8
      a = categorical(softmax(a_logits))

      # return mutant
      seq[:,i] = a
    
    return seq

  def design_semigreedy(self, iters=100, tries=10, dropout=False, use_plddt=True,
                        save_best=True, seq_logits=None, e_tries=None, **kwargs):

    '''semigreedy search'''    
    if e_tries is None: e_tries = tries

    plddt = None
    seq = self._params["seq"].argmax(-1)
    model_flags = {k:kwargs.pop(k,None) for k in ["num_models","sample_models","models"]}
    verbose = kwargs.pop("verbose",1)

    # optimize!
    if verbose: print("Running semigreedy optimization...")
    for i in range(iters):

      buff = []
      model_nums = self._get_model_nums(**model_flags)
      num_tries = (tries+(e_tries-tries)*((i+1)/iters))
      for t in range(int(num_tries)):
        mut_seq = self._mutate(seq, plddt, logits=seq_logits)
        aux = self.predict(mut_seq, return_aux=True, model_nums=model_nums,
                           verbose=False, **kwargs)
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

  def design_pssm_semigreedy(self, soft_iters=300, hard_iters=32, tries=10, e_tries=None,
                             ramp_recycles=True, ramp_models=True, **kwargs):

    verbose = kwargs.get("verbose",1)

    # stage 1: logits -> softmax(logits)
    if soft_iters > 0:
      self.design_3stage(soft_iters, 0, 0, ramp_recycles=ramp_recycles, **kwargs)
      self._tmp["seq_logits"] = kwargs["seq_logits"] = self.aux["seq"]["logits"]

    # stage 2: semi_greedy
    if hard_iters > 0:
      kwargs["dropout"] = False
      if ramp_models:
        num_models = len(kwargs.get("models",self._model_names))
        iters = hard_iters
        for m in range(num_models):
          if verbose and m > 0: print(f'Increasing number of models to {m+1}.')

          kwargs["num_models"] = m + 1
          kwargs["save_best"] = (m + 1) == num_models
          self.design_semigreedy(iters, tries=tries, e_tries=e_tries, **kwargs)
          if m < 2: iters = iters // 2
      else:
        self.design_semigreedy(hard_iters, tries=tries, e_tries=e_tries, **kwargs)


  # ---------------------------------------------------------------------------------
  # experimental optimizers (not extensively evaluated)
  # ---------------------------------------------------------------------------------

  def _design_mcmc(self, steps=1000, half_life=200, T_init=0.01, mutation_rate=1,
                   use_plddt=True, seq_logits=None, save_best=True, **kwargs):
    '''
    MCMC with simulated annealing
    ----------------------------------------
    steps = number for steps for the MCMC trajectory
    half_life = half-life for the temperature decay during simulated annealing
    T_init = starting temperature for simulated annealing. Temperature is decayed exponentially
    mutation_rate = number of mutations at each MCMC step
    '''

    # code borrowed from: github.com/bwicky/oligomer_hallucination

    # gather settings
    verbose = kwargs.pop("verbose",1)
    model_flags = {k:kwargs.pop(k,None) for k in ["num_models","sample_models","models"]}

    # initalize
    plddt, best_loss, current_loss = None, np.inf, np.inf 
    current_seq = self._params["seq"].argmax(-1)

    # run!
    if verbose: print("Running MCMC with simulated annealing...")
    for i in range(steps):

      # update temperature
      T = T_init * (np.exp(np.log(0.5) / half_life) ** i) 

      # mutate sequence
      if i == 0: mut_seq = current_seq
      else:      mut_seq = self._mutate(current_seq, plddt, seq_logits, mutation_rate)

      # get loss
      model_nums = self._get_model_nums(**model_flags)
      aux = self.predict(mut_seq, return_aux=True, verbose=False, model_nums=model_nums, **kwargs)
      loss = aux["log"]["loss"]
  
      # decide
      delta = loss - current_loss
      if i == 0 or delta < 0 or np.random.uniform() < np.exp( -delta / T):

        # accept
        (current_seq,current_loss) = (mut_seq,loss)
        
        if use_plddt:
          plddt = aux["all"]["plddt"].mean(0)
          plddt = plddt[self._target_len:] if self.protocol == "binder" else plddt[:self._len]
        
        if loss < best_loss:
          (best_loss, self._k) = (loss, i)
          self.set_seq(seq=current_seq)
          self._save_results(save_best=save_best, verbose=verbose)