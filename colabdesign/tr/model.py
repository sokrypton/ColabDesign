import random, os
import numpy as np
import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt

from colabdesign.shared.utils import update_dict, Key, dict_to_str
from colabdesign.shared.protein import _np_get_6D_binned
from colabdesign.tr.trrosetta import TrRosetta, get_model_params

# borrow some stuff from AfDesign
from colabdesign.af.prep import prep_pdb, prep_pos
from colabdesign.af.alphafold.common import protein, residue_constants
ORDER_RESTYPE = {v: k for k, v in residue_constants.restype_order.items()}

try:
  from jax.example_libraries.optimizers import sgd, adam
except:
  from jax.experimental.optimizers import sgd, adam

class mk_trdesign_model():
  def __init__(self, protocol="fixbb", num_models=1,
               sample_models=True, data_dir="params/tr",
               loss_callback=None):
    
    assert protocol in ["fixbb","hallucination","partial"]

    self.protocol = protocol
    self._data_dir = "." if os.path.isfile("models/model_xaa.npy") else data_dir
    self._loss_callback = loss_callback

    # set default options
    self._opt = {"temp":1.0, "soft":1.0, "hard":1.0,
                 "models":num_models,"sample_models":sample_models,
                 "weights":{}, "lr":1.0, "bias":0.0, "alpha":1.0}
                
    self.params = {"seq":None}

    # setup model
    self._runner = TrRosetta()
    self._model_params = [get_model_params(f"{self._data_dir}/models/model_xa{k}.npy") for k in list("abcde")]
    self._grad_fn, self._fn = [jax.jit(x) for x in self._get_model()]

    if protocol in ["hallucination","partial"]:
      self.bkg_model = TrRosetta(bkg_model=True)
  
  def _get_model(self):

    def _get_seq(params, opt):
      seq = {"input":params["seq"]}

      # straight-through/reparameterization
      seq["logits"] = seq["input"] * opt["alpha"] + opt["bias"]
      seq["soft"] = jax.nn.softmax(seq["logits"] / opt["temp"])
      seq["hard"] = jax.nn.one_hot(seq["soft"].argmax(-1), 20)
      seq["hard"] = jax.lax.stop_gradient(seq["hard"] - seq["soft"]) + seq["soft"]

      # create pseudo sequence
      seq["pseudo"] = opt["soft"] * seq["soft"] + (1-opt["soft"]) * seq["input"]
      seq["pseudo"] = opt["hard"] * seq["hard"] + (1-opt["hard"]) * seq["pseudo"]

      if self.protocol in ["partial"] and "pos" in opt:
        pos = opt["pos"]
        seq_ref = jax.nn.one_hot(self._batch["aatype"],20)
        seq = jax.tree_map(lambda x: jnp.where(opt["fix_seq"],x.at[pos,:].set(seq_ref),x), seq)
      return seq
    
    def _get_loss(outputs, opt):
      aux = {"outputs":outputs, "losses":{}}
      log_p = jax.tree_map(jax.nn.log_softmax, outputs)

      # bkg loss
      if self.protocol in ["hallucination","partial"]:
        p = jax.tree_map(jax.nn.softmax, outputs)
        log_q = jax.tree_map(jax.nn.log_softmax, self._bkg_feats)
        aux["losses"]["bkg"] = {}
        for k in ["dist","omega","theta","phi"]:
          aux["losses"]["bkg"][k] = -(p[k]*(log_p[k]-log_q[k])).sum(-1).mean()

      # cce loss
      if self.protocol in ["fixbb","partial"]:
        if self.protocol in ["partial"] and "pos" in opt:
          pos = opt["pos"]
          log_p = jax.tree_map(lambda x:x[:,pos][pos,:], log_p)

        q = self._feats
        aux["losses"]["cce"] = {}
        for k in ["dist","omega","theta","phi"]:
          aux["losses"]["cce"][k] = -(q[k]*log_p[k]).sum(-1).mean()

      if self._loss_callback is not None:
        aux["losses"].update(self._loss_callback(outputs))

      # weighted loss
      w = opt["weights"]
      tree_multi = lambda x,y: jax.tree_map(lambda a,b:a*b, x,y)
      losses = {k:(tree_multi(v,w[k]) if k in w else v) for k,v in aux["losses"].items()}
      loss = sum(jax.tree_leaves(losses))
      return loss, aux

    def _model(params, model_params, opt):
      seq = _get_seq(params, opt)
      outputs = self._runner(seq["pseudo"], model_params)
      loss, aux = _get_loss(outputs, opt)
      aux.update({"seq":seq,"opt":opt})
      return loss, aux

    return jax.value_and_grad(_model, has_aux=True, argnums=0), _model
  
  def prep_inputs(self, pdb_filename=None, chain=None, length=None,
                  pos=None, fix_seq=False, **kwargs):
    '''
    prep inputs for TrDesign
    '''    
    if self.protocol in ["fixbb", "partial"]:
      # parse PDB file and return features compatible with TrRosetta
      pdb = prep_pdb(pdb_filename, chain, for_alphafold=False)
      self._batch = pdb["batch"]

      if self.protocol in ["partial"] and pos is not None:
        p = prep_pos(pos, **pdb["idx"])
        self._batch = jax.tree_map(lambda x:x[p], self._batch)
        self._opt["pos"] = np.arange(len(p))
        self._opt["fix_seq"] = fix_seq

      self._feats = _np_get_6D_binned(self._batch["all_atom_positions"],
                                      self._batch["all_atom_mask"])

      self._len = len(self._batch["aatype"])
      self._opt["weights"]["cce"] = {"dist":1/6,"omega":1/6,"theta":2/6,"phi":2/6}

    if self.protocol in ["hallucination", "partial"]:
      # compute background distribution
      if length is not None: self._len = length
      self._bkg_feats = []
      key = jax.random.PRNGKey(0)
      for n in range(1,6):
        model_params = get_model_params(f"{self._data_dir}/bkgr_models/bkgr0{n}.npy")
        self._bkg_feats.append(self.bkg_model(model_params, key, self._len))
      self._bkg_feats = jax.tree_map(lambda *x:jnp.stack(x).mean(0), *self._bkg_feats)
      self._opt["weights"]["bkg"] = {"dist":1/6,"omega":1/6,"theta":2/6,"phi":2/6}

    self.restart(**kwargs)

  def set_opt(self, *args, **kwargs):
    update_dict(self.opt, *args, **kwargs)

  def set_weights(self, *args, **kwargs):
    update_dict(self.opt["weights"], *args, **kwargs)
  
  def _init_seq(self, seq=None):
    if seq is None:
      seq = 0.01 * jax.random.normal(self.key(), (self._len,20))
    else:
      if isinstance(seq, str):
        seq = np.array([residue_constants.restype_order.get(aa,-1) for aa in seq])
        seq = np.eye(20)[seq]
      assert seq.shape[-2] == self._len
    update_dict(self.params, {"seq":seq})

  def _init_optimizer(self):
    '''setup which optimizer to use'''      
    self._init_fun, self._update_fun, self._get_params = sgd(1.0)
    self._k = 0
    self._state = self._init_fun(self.params)

  def restart(self, seed=None, opt=None, weights=None, seq=None):
    self.key = Key(seed=seed).get
    self.opt = jax.tree_map(lambda x:x, self._opt) # copy
    self.set_opt(opt)
    self.set_weights(weights)
    self._init_seq(seq)
    self._init_optimizer()

    self._best_loss = np.inf
    self._best_aux = None
    
  def run(self, backprop=True):
    '''run model to get outputs, losses and gradients'''
    
    # decide which model params to use
    m = self.opt["models"]
    ns = jnp.arange(5)
    if self.opt["sample_models"] and m != len(ns):
      model_num = jax.random.choice(self.key(),ns,(m,),replace=False)
    else:
      model_num = ns[:m]
    model_num = np.array(model_num).tolist()

    # run in serial
    _loss, _aux, _grad = [],[],[]
    for n in model_num:
      model_params = self._model_params[n]
      if backprop:
        (l,a),g = self._grad_fn(self.params, model_params, self.opt)
        _grad.append(g)
      else:
        l,a = self._fn(self.params, model_params, self.opt)
      _loss.append(l)
      _aux.append(a)
    
    # average results
    if len(model_num) > 1:
      _loss = jnp.asarray(_loss).mean()
      _aux = jax.tree_map(lambda *v: jnp.stack(v).mean(0), *_aux)
      if backprop: _grad = jax.tree_map(lambda *v: jnp.stack(v).mean(0), *_grad)
    else:
      _loss,_aux = _loss[0],_aux[0]
      if backprop: _grad = _grad[0] 

    if not backprop:
      _grad = jax.tree_map(jnp.zeros_like, self.params)

    # update
    self.loss = _loss
    self.aux = _aux
    self.aux["model_num"] = model_num
    self.grad = _grad

  def step(self, backprop=True, callback=None,
           save_best=True, verbose=1):
    self.run(backprop=backprop)
    if callback is not None: callback(self)

    # normalize gradient
    g = self.grad["seq"]
    gn = jnp.linalg.norm(g,axis=(-1,-2),keepdims=True)
    self.grad["seq"] *= jnp.sqrt(self._len)/(gn+1e-7)

    # set learning rate
    self.grad = jax.tree_map(lambda x:x*self.opt["lr"], self.grad)

    # apply gradient
    self._state = self._update_fun(self._k, self.grad, self._state)
    self.params = self._get_params(self._state)    

    # increment
    self._k += 1

    # save results
    if save_best and self.loss < self._best_loss:
      self._best_loss = self.loss
      self._best_aux = self.aux

    # print
    if verbose and (self._k % verbose) == 0:
      x = self.get_loss(get_best=False)
      x["models"] = self.aux["model_num"]
      print(dict_to_str(x, print_str=f"{self._k}", keys=["models"]))

  def design(self, iters=100, opt=None, weights=None, save_best=True, verbose=1):
    self.set_opt(opt)
    self.set_weights(weights)
    for _ in range(iters):
      self.step(save_best=save_best, verbose=verbose)

  def predict(self):
    '''make prediction'''
    self.run(backprop=False)
  
  def plot(self, mode="preds", get_best=True, dpi=100):
    '''plot predictions'''

    assert mode in ["preds","feats","bkg_feats"]
    if mode == "preds":
      aux = self.aux if (self._best_aux is None or not get_best) else self._best_aux
      x = aux["outputs"]
    elif mode == "feats":
      x = self._feats
    elif mode == "bkg_feats":
      x = self._bkg_feats

    x = jax.tree_map(np.asarray, x)

    plt.figure(figsize=(4*4,4), dpi=dpi)
    for n,k in enumerate(["theta","phi","dist","omega"]):
      v = x[k]
      plt.subplot(1,4,n+1)
      plt.title(k)
      plt.imshow(v.argmax(-1),cmap="binary")
    plt.show()
    
  def get_loss(self, k=None, get_best=True):
    aux = self.aux if (self._best_aux is None or not get_best) else self._best_aux
    if k is None:
      return {k:self.get_loss(k,get_best=get_best) for k in aux["losses"].keys()}
    losses = aux["losses"][k]
    weights = aux["opt"]["weights"][k]
    weighted_losses = jax.tree_map(lambda l,w:l*w, losses, weights)
    return float(sum(jax.tree_leaves(weighted_losses)))

  def get_seq(self, get_best=True):
    aux = self.aux if (self._best_aux is None or not get_best) else self._best_aux
    x = np.array(aux["seq"]["pseudo"].argmax(-1))
    return "".join([ORDER_RESTYPE[a] for a in x])
    
  def af_callback(self, weight=1.0, seed=None):
    
    def callback(af_model):
      
      # copy [opt]ions from afdesign
      for k in ["soft","temp","hard","pos","fix_seq","bias","models"]:
        if k in self.opt and k in af_model.opt:
          self.opt[k] = af_model.opt[k]

      # update sequence input
      self.params["seq"] = af_model.params["seq"][0]
      
      # run trdesign
      self.run(backprop = weight > 0)
      
      # add gradients
      af_model.grad["seq"] += weight * self.grad["seq"]
      
      # add loss
      af_model.loss += weight * self.loss
        
      # for verbose printout
      if self.protocol in ["hallucination","partial"]:
        af_model.aux["losses"]["TrD_bkg"] = self.get_loss("bkg")
      if self.protocol in ["fixbb","partial"]:
        af_model.aux["losses"]["TrD_cce"] = self.get_loss("cce")
      
    self.restart(seed=seed)
    return callback