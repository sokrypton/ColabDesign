import random
import numpy as np
import jax
import jax.numpy as jnp

from tr.src.utils import update_dict
from tr.src.trrosetta import TrRosetta, get_model_params

# borrow some stuff from AfDesign
from af.src.misc import _np_get_6D
from af.src.prep import prep_pdb, prep_pos
from alphafold.common import protein, residue_constants
ORDER_RESTYPE = {v: k for k, v in residue_constants.restype_order.items()}

class mk_trdesign_model():
  def __init__(self, protocol="fixbb", model_num=1, model_sample=True,
               seed=None, data_dir="."):
    
    self.protocol = protocol

    # set default options
    self.opt = {"temp":1.0, "soft":1.0, "hard":1.0,
                "weights":{"dist":0.25,"omega":0.25,
                           "theta":0.25,"phi":0.25},
                "model":{"num":model_num,"sample":model_sample}}
                
    self.params = {"seq":None}
    self._seed = random.randint(0,2147483647) if seed is None else seed
    self._key = jax.random.PRNGKey(self._seed)

    # setup model
    self._runner = TrRosetta()
    self._model_params = [get_model_params(f"{data_dir}/models/model_xa{k}.npy") for k in list("abcde")]
    self._grad_fn, self._fn = [jax.jit(x) for x in self._get_model()]
  
  def _get_model(self):

    def _get_seq(params, opt):
      seq = {"input":params["seq"]}

      # straight-through/reparameterization
      seq["logits"] = 2.0 * seq["input"]
      seq["soft"] = jax.nn.softmax(seq["logits"] / opt["temp"])
      seq["hard"] = jax.nn.one_hot(seq["soft"].argmax(-1), 20)
      seq["hard"] = jax.lax.stop_gradient(seq["hard"] - seq["soft"]) + seq["soft"]

      # create pseudo sequence
      seq["pseudo"] = opt["soft"] * seq["soft"] + (1-opt["soft"]) * seq["input"]
      seq["pseudo"] = opt["hard"] * seq["hard"] + (1-opt["hard"]) * seq["pseudo"]

      if self.protocol in ["partial"]:
        p = opt["pos"]
        seq_ref = jax.nn.one_hot(self._wt_aatype,20)
        seq = jax.tree_map(lambda x:x.at[:,p,:].set(seq_ref), seq)

      return seq
    
    def _get_loss(outputs, opt):
      aux = {"outputs":outputs, "losses":{}}
      if self.protocol in ["fixbb"]:
        log_p = jax.tree_map(jax.nn.log_softmax, outputs)
        for k in ["dist","omega","theta","phi"]:
          q = self.feats[k]
          aux["losses"][k] = -(q*log_p[k]).sum(-1).mean()

      # weighted loss
      loss = []
      losses_keys = list(aux["losses"].keys())
      for k in losses_keys:
        if k in opt["weights"]:
          loss.append(aux["losses"][k] * opt["weights"][k])
        else:
          aux["losses"].pop(k)
      return sum(loss), aux

    def _model(params, model_params, opt):
      seq = _get_seq(params, opt)
      outputs = self._runner(seq["pseudo"], model_params)

      loss, aux = _get_loss(outputs, opt)
      return loss, aux

    return jax.value_and_grad(_model, has_aux=True, argnums=0), _model
  
  def prep_inputs(self, pdb_filename=None, chain=None, length=None, pos=None):
    '''Parse PDB file and return features compatible with TrRosetta'''
    if pdb_filename is not None:
      pdb = prep_pdb(pdb_filename, chain, for_alphafold=False)
      ref = _np_get_6D(pdb["batch"]["all_atom_positions"],
                       pdb["batch"]["all_atom_mask"],
                       use_jax=False, for_trrosetta=True)
      ref = jax.tree_map(jnp.squeeze,ref)

      def mtx2bins(x_ref, start, end, nbins, mask):
        bins = np.linspace(start, end, nbins)
        x_true = np.digitize(x_ref, bins).astype(np.uint8)
        x_true[mask] = 0
        return np.eye(nbins+1)[x_true][...,:-1]

      mask = (ref["dist"] > 20) | (np.eye(ref["dist"].shape[0]) == 1)
      self.feats = {"dist":mtx2bins(ref["dist"],     2.0,  20.0,  37, mask=mask),
                    "omega":mtx2bins(ref["omega"], -np.pi, np.pi, 25,  mask=mask),
                    "theta":mtx2bins(ref["theta"], -np.pi, np.pi, 25,  mask=mask),
                    "phi":  mtx2bins(ref["phi"],      0.0, np.pi, 13,  mask=mask)}
      
      self._wt_seq = pdb["batch"]["aatype"]
      self._len = len(self._wt_seq)      

      if self.params["seq"] is None:
        self.params["seq"] = np.zeros((self._len,20))
    
  def run(self, seq=None, params=None, opt=None, weights=None, backprop=True):
    # update settings if defined
    update_dict(self.params, {"seq":seq})
    update_dict(self.params, params)
    update_dict(self.opt, opt)
    update_dict(self.opt, {"weights":weights})

    # decide which model params to use
    m = self.opt["model"]["num"]
    ns = jnp.arange(5)
    if self.opt["model"]["sample"] and m != len(ns):
      self._key, key = jax.random.split(self._key)
      model_num = jax.random.choice(key,ns,(m,),replace=False)
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
      _aux = jax.tree_map(lambda *v: jnp.asarray(v).mean(0), *_aux)
      if backprop:
        _grad = jax.tree_map(lambda *v: jnp.asarray(v).mean(0), *_grad)
      else:
        _grad = None
    else:
      _loss,_aux = _loss[0],_aux[0]
      _grad = _grad[0] if backprop else None
    
    return {"loss":_loss,"grad":_grad,
            "aux":_aux,"model_num":model_num}
