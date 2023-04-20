import random, os
import jax
import jax.numpy as jnp
import numpy as np
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.shared.utils import copy_dict, update_dict, Key, dict_to_str, to_float, softmax, categorical, to_list, copy_missing

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
              seq=None, mode=None, keep_history=False, reset_opt=True, **kwargs):   
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
      copy_missing(self.opt, self._opt)
      self.opt = copy_dict(self._opt)
      if hasattr(self,"aux"): del self.aux
    
    if not keep_history:
      # initialize trajectory
      self._tmp = {"traj":{"seq":[],"xyz":[],"plddt":[],"pae":[]},
                   "log":[],"best":{}}

    # update options/settings (if defined)
    self.set_opt(opt)
    self.set_weights(weights)
  
    # initialize sequence
    self.set_seed(seed)
    self.set_seq(seq=seq, mode=mode, **kwargs)

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
      model_nums = np.random.choice(ns,(m,),replace=False)
    else:
      model_nums = ns[:m]
    return model_nums   

  def run(self, num_recycles=None, num_models=None, sample_models=None, models=None,
          backprop=True, callback=None, model_nums=None, return_aux=False):
    '''run model to get outputs, losses and gradients'''
    
    # pre-design callbacks
    for fn in self._callbacks["design"]["pre"]: fn(self)

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

    # update aux (average outputs)
    def avg_or_first(x):
      if np.issubdtype(x.dtype, np.integer): return x[0]
      else: return x.mean(0)

    self.aux = jax.tree_map(avg_or_first, auxs)
    self.aux["atom_positions"] = auxs["atom_positions"][0]
    self.aux["all"] = auxs
    
    # post-design callbacks
    for fn in (self._callbacks["design"]["post"] + to_list(callback)): fn(self)

    # update log
    self.aux["log"] = {**self.aux["losses"]}
    self.aux["log"]["plddt"] = 1 - self.aux["log"]["plddt"]
    for k in ["loss","i_ptm","ptm"]: self.aux["log"][k] = self.aux[k]
    for k in ["hard","soft","temp"]: self.aux["log"][k] = self.opt[k]

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
    self.aux["log"].update({"recycles":int(self.aux["num_recycles"]),
                            "models":model_nums})
    
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
    a = self._args
    mode = a["recycle_mode"]
    if num_recycles is None:
      num_recycles = self.opt["num_recycles"]

    if mode in ["backprop","add_prev"]:
      # recycles compiled into model, only need single-pass
      aux = self._single(model_params, backprop)
    
    else:
      L = self._inputs["residue_index"].shape[0]
      
      # intialize previous
      if "prev" not in self._inputs or a["clear_prev"]:
        prev = {'prev_msa_first_row': np.zeros([L,256]),
                'prev_pair': np.zeros([L,L,128])}

        if a["use_initial_guess"] and "batch" in self._inputs:
          prev["prev_pos"] = self._inputs["batch"]["all_atom_positions"] 
        else:
          prev["prev_pos"] = np.zeros([L,37,3])

        if a["use_dgram"]:
          # TODO: add support for initial_guess + use_dgram
          prev["prev_dgram"] = np.zeros([L,L,64])

        if a["use_initial_atom_pos"]:
          if "batch" in self._inputs:
            self._inputs["initial_atom_pos"] = self._inputs["batch"]["all_atom_positions"] 
          else:
            self._inputs["initial_atom_pos"] = np.zeros([L,37,3])              
      
      self._inputs["prev"] = prev
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
        if a["use_initial_atom_pos"]:
          self._inputs["initial_atom_pos"] = aux["prev"]["prev_pos"]                

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

    # save results
    self._save_results(save_best=save_best, verbose=verbose)

    # increment
    self._k += 1

  def _print_log(self, print_str=None, aux=None):
    if aux is None: aux = self.aux
    keys = ["models","recycles","hard","soft","temp","seqid","loss",
            "seq_ent","mlm","helix","pae","i_pae","exp_res","con","i_con",
            "sc_fape","sc_rmsd","dgram_cce","fape","plddt","ptm"]
    
    if "i_ptm" in aux["log"]:
      if len(self._lengths) > 1:
        keys.append("i_ptm")
      else:
        aux["log"].pop("i_ptm")

    print(dict_to_str(aux["log"], filt=self.opt["weights"],
                      print_str=print_str, keys=keys+["rmsd"], ok=["plddt","rmsd"]))

  def _save_results(self, aux=None, save_best=False,
                    best_metric=None, metric_higher_better=False,
                    verbose=True):
    if aux is None: aux = self.aux    
    self._tmp["log"].append(aux["log"])    
    if (self._k % self._args["traj_iter"]) == 0:
      # update traj
      traj = {"seq":   aux["seq"]["pseudo"],
              "xyz":   aux["atom_positions"][:,1,:],
              "plddt": aux["plddt"],
              "pae":   aux["pae"]}
      for k,v in traj.items():
        if len(self._tmp["traj"][k]) == self._args["traj_max"]:
          self._tmp["traj"][k].pop(0)
        self._tmp["traj"][k].append(v)

    # save best
    if save_best:
      if best_metric is None:
        best_metric = self._args["best_metric"]
      metric = float(aux["log"][best_metric])
      if self._args["best_metric"] in ["plddt","ptm","i_ptm","seqid","composite"] or metric_higher_better:
        metric = -metric
      if "metric" not in self._tmp["best"] or metric < self._tmp["best"]["metric"]:
        self._tmp["best"]["aux"] = copy_dict(aux)
        self._tmp["best"]["metric"] = metric

    if verbose and ((self._k+1) % verbose) == 0:
      self._print_log(f"{self._k+1}", aux=aux)

  def predict(self, seq=None, bias=None,
              num_models=None, num_recycles=None, models=None, sample_models=False,
              dropout=False, hard=True, soft=False, temp=1,
              return_aux=False, verbose=True,  seed=None, **kwargs):
    '''predict structure for input sequence (if provided)'''

    def load_settings():    
      if "save" in self._tmp:
        [self.opt, self._args, self._params, self._inputs] = self._tmp.pop("save")

    def save_settings():
      load_settings()
      self._tmp["save"] = [copy_dict(x) for x in [self.opt, self._args, self._params, self._inputs]]

    save_settings()

    # set seed if defined
    if seed is not None: self.set_seed(seed)

    # set [seq]uence/[opt]ions
    if seq is not None: self.set_seq(seq=seq, bias=bias)    
    self.set_opt(hard=hard, soft=soft, temp=temp, dropout=dropout, pssm_hard=True)
    self.set_args(shuffle_first=False)
    
    # run
    self.run(num_recycles=num_recycles, num_models=num_models,
             sample_models=sample_models, models=models, backprop=False, **kwargs)
    if verbose: self._print_log("predict")

    load_settings()

    # return (or save) results
    if return_aux: return self.aux

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
    i_prob[np.isnan(i_prob)] = 0
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
      a_logits = logits - np.eye(self._args["alphabet_size"])[seq[:,i]] * 1e8
      a = categorical(softmax(a_logits))

      # return mutant
      seq[:,i] = a
    
    return seq

  def design_semigreedy(self, iters=100, tries=10, dropout=False,
                        save_best=True, seq_logits=None, e_tries=None, **kwargs):

    '''semigreedy search'''    
    if e_tries is None: e_tries = tries

    # get starting sequence
    if hasattr(self,"aux"):
      seq = self.aux["seq"]["logits"].argmax(-1)
    else:
      seq = (self._params["seq"] + self._inputs["bias"]).argmax(-1)

    # bias sampling towards the defined bias
    if seq_logits is None: seq_logits = 0
    
    model_flags = {k:kwargs.pop(k,None) for k in ["num_models","sample_models","models"]}
    verbose = kwargs.pop("verbose",1)

    # get current plddt
    aux = self.predict(seq, return_aux=True, verbose=False, **model_flags, **kwargs)
    plddt = self.aux["plddt"]
    plddt = plddt[self._target_len:] if self.protocol == "binder" else plddt[:self._len]

    # optimize!
    if verbose:
      print("Running semigreedy optimization...")
    
    for i in range(iters):
      buff = []
      model_nums = self._get_model_nums(**model_flags)
      num_tries = (tries+(e_tries-tries)*((i+1)/iters))
      for t in range(int(num_tries)):
        mut_seq = self._mutate(seq=seq, plddt=plddt,
                               logits=seq_logits + self._inputs["bias"])
        aux = self.predict(seq=mut_seq, return_aux=True, model_nums=model_nums, verbose=False, **kwargs)
        buff.append({"aux":aux, "seq":np.array(mut_seq)})

      # accept best
      losses = [x["aux"]["loss"] for x in buff]
      best = buff[np.argmin(losses)]
      self.aux, seq = best["aux"], jnp.array(best["seq"])
      self.set_seq(seq=seq, bias=self._inputs["bias"])
      self._save_results(save_best=save_best, verbose=verbose)

      # update plddt
      plddt = best["aux"]["plddt"]
      plddt = plddt[self._target_len:] if self.protocol == "binder" else plddt[:self._len]
      self._k += 1

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
                   seq_logits=None, save_best=True, **kwargs):
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

    # initialize
    plddt, best_loss, current_loss = None, np.inf, np.inf 
    current_seq = (self._params["seq"] + self._inputs["bias"]).argmax(-1)
    if seq_logits is None: seq_logits = 0

    # run!
    if verbose: print("Running MCMC with simulated annealing...")
    for i in range(steps):

      # update temperature
      T = T_init * (np.exp(np.log(0.5) / half_life) ** i) 

      # mutate sequence
      if i == 0:
        mut_seq = current_seq
      else:
        mut_seq = self._mutate(seq=current_seq, plddt=plddt,
                               logits=seq_logits + self._inputs["bias"],
                               mutation_rate=mutation_rate)

      # get loss
      model_nums = self._get_model_nums(**model_flags)
      aux = self.predict(seq=mut_seq, return_aux=True, verbose=False, model_nums=model_nums, **kwargs)
      loss = aux["log"]["loss"]
  
      # decide
      delta = loss - current_loss
      if i == 0 or delta < 0 or np.random.uniform() < np.exp( -delta / T):

        # accept
        (current_seq,current_loss) = (mut_seq,loss)
        
        plddt = aux["all"]["plddt"].mean(0)
        plddt = plddt[self._target_len:] if self.protocol == "binder" else plddt[:self._len]
        
        if loss < best_loss:
          (best_loss, self._k) = (loss, i)
          self.set_seq(seq=current_seq, bias=self._inputs["bias"])
          self._save_results(save_best=save_best, verbose=verbose)
