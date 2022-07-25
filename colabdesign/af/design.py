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

    # initialize sequence
    self.key = Key(seed=seed).get
    self.set_seq(seq=seq, **kwargs)

    # setup optimizer
    if isinstance(optimizer, str):
      if optimizer == "adam": (optimizer,lr) = (adam,0.02)
      if optimizer == "sgd":  (optimizer,lr) = (sgd,0.1)
      self.opt["lr"] = lr

    self._init_fun, self._update_fun, self._get_params = optimizer(1.0)
    self._state = self._init_fun(self.params)
    self._k = 0

    if not keep_history:
      # initialize trajectory
      self._traj = {"log":[],"seq":[],"xyz":[],"plddt":[],"pae":[]}
      self._best = {"metric":np.inf}

    # set crop length
    if self._args["crop_len"] is None:
      self._args["crop_len"] = self._inputs["residue_index"].shape[-1]
      
  def run(self, model=None, backprop=True, crop=True, callback=None):
    '''run model to get outputs, losses and gradients'''

    # crop inputs 
    (L, max_L) = (self._inputs["residue_index"].shape[-1], self._args["crop_len"])
    if crop and max_L < L:
      crop_mode = self._args["crop_mode"]
    
      if crop_mode == "slide":
        i = jax.random.randint(self.key(),[],0,L-max_L)
        p = np.arange(i,i+max_L)
      
      if crop_mode == "roll":
        i = jax.random.randint(self.key(),[],0,L)
        p = np.sort(np.roll(np.arange(L),L-i)[:max_L])

      # if crop_mode == "dist":
      # TODO: shihao
    
    else:
      p = np.arange(L) 
    
    self.opt["crop_pos"] = p
    
    # decide which model params to use
    ns,ns_name = [],[]
    for n,name in enumerate(self._model_names):
      if "openfold" in name:
        if self._args["use_openfold"]:  ns.append(n); ns_name.append(name)
      elif self._args["use_alphafold"]: ns.append(n); ns_name.append(name)

    # sub select number of model params
    if model is not None:
      model = [model] if isinstance(model,int) else list(model)
      ns = [ns[n if isinstance(n,int) else ns_name.index(n)] for n in model]
    
    ns = jnp.array(ns)
    m = min(self.opt["models"],len(ns))
    if self.opt["sample_models"] and m != len(ns):
      model_num = jax.random.choice(self.key(),ns,(m,),replace=False)
    else:
      model_num = ns[:m]      
    model_num = np.array(model_num).tolist()

    # loop through model params
    aux = []
    for n in model_num:
      p = self._model_params[n]
      aux.append(self._recycle(p, backprop=backprop))
    aux = jax.tree_map(lambda *x: jnp.stack(x), *aux)

    # update aux
    self.aux = jax.tree_map(lambda x:x[0], aux)
    self.aux["all"] = aux

    # average losses and gradients
    self.aux["loss"] = aux["loss"].mean()
    self.aux["losses"] = jax.tree_map(lambda x: x.mean(0), aux["losses"])
    self.aux["grad"] = jax.tree_map(lambda x: x.mean(0), aux["grad"])

    # callback
    if callback is not None: callback(self)

    # update log
    self.aux["log"] = {**self.aux["losses"], "loss":self.aux["loss"]}
    self.aux["log"].update({k:self.opt[k] for k in ["hard","soft","temp"]})

    # compute sequence recovery
    if self.protocol in ["fixbb","partial"] or (self.protocol == "binder" and self._args["redesign"]):
      aatype = self.aux["seq"]["pseudo"].argmax(-1)
      if self.protocol == "partial" and "pos" in self.opt:
        aatype = aatype[...,self.opt["pos"]]
      self.aux["log"]["seqid"] = (aatype == self._wt_aatype).mean()

    self.aux["log"] = to_float(self.aux["log"])
    self.aux["log"].update({"recycles":self.aux["recycles"], "models":model_num})
    
  def _single(self, model_params, backprop=True):
    '''single pass through the model'''
    flags  = [self.params, model_params, self._inputs, self.key(), self.opt]
    if backprop:
      (loss, aux), grad = self._grad_fn(*flags)
    else:
      loss, aux = self._fn(*flags)
      grad = jax.tree_map(jnp.zeros_like, self.params)
    aux.update({"loss":loss,"grad":grad})
    return aux

  def _recycle(self, model_params, backprop=True):   
    '''multiple passes through the model (aka recycle)'''

    mode = self._args["recycle_mode"]
    if mode in ["backprop","add_prev"]:
      # recycles compiled into model, only need single-pass
      recycles = self.opt["recycles"] = self._runner.config.model.num_recycle
      aux = self._single(model_params, backprop)
    
    else:
      # configure number of recycle to run
      recycles = self.opt["recycles"]
      if mode == "average":
        # run recycles manually, average gradients
        if "crop_pos" in self.opt: L = self.opt["crop_pos"].shape[0]
        else: L = self._inputs["residue_index"].shape[-1]
        self._inputs["prev"] = {'prev_msa_first_row': np.zeros([L,256]),
                                'prev_pair': np.zeros([L,L,128]),
                                'prev_pos': np.zeros([L,37,3])}
        grad = []
        for _ in range(recycles+1):
          aux = self._single(model_params, backprop)
          grad.append(aux["grad"])
          self._inputs["prev"] = aux["prev"]
        aux["grad"] = jax.tree_map(lambda *x: jnp.stack(x).mean(0), *grad)
      
      elif mode == "sample":
        # randomly select number of recycles to run
        self.set_opt(recycles=jax.random.randint(self.key(),[],0,recycles+1))
        aux = self._single(model_params, backprop)
        (self.opt["recycles"],recycles) = (recycles,self.opt["recycles"])
      
      else:
        aux = self._single(model_params, backprop)
    
    aux["recycles"] = recycles
    return aux

  def step(self, lr_scale=1.0, model=None, backprop=True, crop=True,
           callback=None, save_best=False, verbose=1):
    '''do one step of gradient descent'''
    
    # run
    self.run(model=model, backprop=backprop, crop=crop, callback=callback)

    # normalize gradient
    g = self.aux["grad"]["seq"]
    gn = jnp.linalg.norm(g,axis=(-1,-2),keepdims=True)
    
    eff_len = (jnp.square(g).sum(-1,keepdims=True) > 0).sum(-2,keepdims=True)
    self.aux["grad"]["seq"] *= jnp.sqrt(eff_len)/(gn+1e-7)

    # set learning rate
    lr = self.opt["lr"] * lr_scale
    self.aux["grad"] = jax.tree_map(lambda x:x*lr, self.aux["grad"])

    # apply gradient
    self._state = self._update_fun(self._k, self.aux["grad"], self._state)
    self.params = self._get_params(self._state)

    # increment
    self._k += 1

    # save results
    self.save_results(save_best=save_best, verbose=verbose)

  def save_results(self, save_best=False, verbose=1):
    
    # save trajectory
    traj = {"seq":   self.aux["seq"]["pseudo"],
            "xyz":   self.aux["atom_positions"][:,1,:],
            "plddt": self.aux["plddt"],
            "pae":   self.aux["pae"]}
    traj = jax.tree_map(np.array, traj)
    traj["log"] = self.aux["log"]
    
    for k,v in traj.items(): self._traj[k].append(v)

    # save best result
    if save_best:
      metric = self.aux["log"][self._args["best_metric"]]
      if metric < self._best["metric"]:
        self._best.update({"metric":metric, "aux":self.aux})

    if verbose and (self._k % verbose) == 0:
      # preferred order
      keys = ["models","recycles","hard","soft","temp","seqid","loss",
              "msa_ent","plddt","pae","helix","con","i_pae","i_con",
              "sc_fape","sc_rmsd","dgram_cce","fape","rmsd"]

      print(dict_to_str(self.aux["log"], filt=self.opt["weights"],
                        print_str=f"{self._k}", keys=keys, ok="rmsd"))

  # ---------------------------------------------------------------------------------
  # example design functions
  # ---------------------------------------------------------------------------------
  def design(self, iters=100,
             soft=None, e_soft=None,
             temp=None, e_temp=None,
             hard=None, e_hard=None,
             opt=None, weights=None, dropout=None, crop=True,
             backprop=True, callback=None, save_best=False, verbose=1):
      
    # update options/settings (if defined)
    self.set_opt(opt, dropout=dropout)
    self.set_weights(weights)

    if soft is None: soft = self.opt["soft"]
    if temp is None: temp = self.opt["temp"]
    if hard is None: hard = self.opt["hard"]
    if e_soft is None: e_soft = soft
    if e_temp is None: e_temp = temp
    if e_hard is None: e_hard = hard
    
    for i in range(iters):
      self.set_opt(soft=(soft+(e_soft-soft)*((i+1)/iters)),
                   hard=(hard+(e_hard-hard)*((i+1)/iters)),
                   temp=(e_temp+(temp-e_temp)*(1-(i+1)/iters)**2))
      
      # decay learning rate based on temperature
      lr_scale = (1 - self.opt["soft"]) + (self.opt["soft"] * self.opt["temp"])
      
      self.step(lr_scale=lr_scale, backprop=backprop, crop=crop,
                callback=callback, save_best=save_best, verbose=verbose)

  def design_logits(self, iters=100, **kwargs):
    '''optimize logits'''
    self.design(iters, soft=0, hard=0, **kwargs)

  def design_soft(self, iters=100, **kwargs):
    ''' optimize softmax(logits/temp)'''
    self.design(iters, soft=1, hard=0, **kwargs)
  
  def design_hard(self, iters=100, **kwargs):
    ''' optimize argmax(logits)'''
    self.design(iters, soft=1, hard=1, **kwargs)

  # ---------------------------------------------------------------------------------
  # experimental
  # ---------------------------------------------------------------------------------

  def design_2stage(self, soft_iters=100, temp_iters=100, hard_iters=50,
                    models=1, dropout=True, **kwargs):
    '''two stage design (soft→hard)'''
    self.set_opt(models=models, sample_models=True) # sample models
    self.design(soft_iters, soft=1, temp=1, dropout=dropout, **kwargs)
    self.design(temp_iters, soft=1, temp=1, dropout=dropout,  e_temp=1e-2, **kwargs)
    self.set_opt(models=len(self._model_params)) # use all models
    self.design(hard_iters, soft=1, temp=1e-2, dropout=False, hard=True, save_best=True, **kwargs)

  def design_3stage(self, soft_iters=300, temp_iters=100, hard_iters=10,
                    models=1, dropout=True, **kwargs):
    '''three stage design (logits→soft→hard)'''
    self.set_opt(models=models, sample_models=True) # sample models
    self.design(soft_iters, soft=0, temp=1,    hard=0, e_soft=1,    dropout=dropout, **kwargs)
    self.design(temp_iters, soft=1, temp=1,    hard=0, e_temp=1e-2, dropout=dropout, **kwargs)
    self.set_opt(models=len(self._model_params)) # use all models
    self.design(hard_iters, soft=1, temp=1e-2, hard=1, dropout=False, save_best=True, **kwargs)

  def design_semigreedy(self, iters=100, tries=20, models=1,
                        use_plddt=True, save_best=True,
                        verbose=1):
    '''semigreedy search'''    
    self.set_opt(hard=True, dropout=False,
                 models=models, sample_models=False)
    if self._k == 0:
      self.run(backprop=False, crop=False)

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
        plddt = np.asarray(self.aux["all"]["plddt"].mean(0))
        if self.protocol == "binder":
          plddt = plddt[self._target_len:]
        else:
          plddt = plddt[:self._len]
      
      buff = []
      for _ in range(tries):
        self.params["seq"] = mut(seq, plddt)
        self.run(backprop=False, crop=False)
        buff.append({"aux":self.aux, "seq":self.params["seq"]})
      
      # accept best      
      buff = buff[np.argmin([x["aux"]["loss"] for x in buff])]
      self.params["seq"], self.aux = buff["seq"], buff["aux"]
      self._k += 1
      self.save_results(save_best=save_best, verbose=verbose)

  # ---------------------------------------------------------------------------------

  def template_predesign(self, iters=100, soft=True, hard=False, 
                         temp=1.0, dropout=True, **kwargs):
    '''use template for predesign stage, then increase template dropout until gone'''
    self.set_opt(hard=hard, soft=soft, dropout=dropout, temp=temp)
    for i in range(iters):
      self.opt["template"]["dropout"] = (i+1)/iters
      self.step(**kwargs)