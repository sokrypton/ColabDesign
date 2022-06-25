import random

import numpy as np
import jax
import jax.numpy as jnp
from tr.src.trrosetta import TrRosetta, get_model_params

class mk_design_model():
  def __init__(self, protocol="fixbb", model_num=1, model_sample=True,
               seed=None, data_dir="."):
    
    self.protocol = protocol

    # set default options
    self.opt = {"temp":1.0, "soft":1.0, "hard":1.0,
                "weights":{},
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
      loss = jnp.square(outputs["dist"]).sum()
      return loss, aux

    def _model(params, model_params, opt):
      seq = _get_seq(params, opt)
      outputs = self._runner(seq["pseudo"], model_params)

      loss, aux = _get_loss(outputs, opt)
      return loss, aux

    return jax.value_and_grad(_model, has_aux=True, argnums=0), _model
  
  def prep_inputs(self, pdb=None, chain=None, length=None, pos=None):
    '''Parse PDB file and return features compatible with TrRosetta'''
    if pdb is not None:

      
      ncac, seq = parse_PDB(pdb,["N","CA","C"], chain=chain)

      # mask gap regions
      if mask_gaps:
        mask = seq != 20
        ncac, seq = ncac[mask], seq[mask]

      N,CA,C = ncac[:,0], ncac[:,1], ncac[:,2]
      CB = extend(C, N, CA, 1.522, 1.927, -2.143)

      dist_ref  = to_len(CB[:,None], CB[None,:])
      omega_ref = to_dih(CA[:,None], CB[:,None], CB[None,:], CA[None,:])
      theta_ref = to_dih( N[:,None], CA[:,None], CB[:,None], CB[None,:])
      phi_ref   = to_ang(CA[:,None], CB[:,None], CB[None,:])

      def mtx2bins(x_ref, start, end, nbins, mask):
        bins = np.linspace(start, end, nbins)
        x_true = np.digitize(x_ref, bins).astype(np.uint8)
        x_true[mask] = 0
        return np.eye(nbins+1)[x_true][...,:-1]

      mask = (dist_ref > 20)
      feats = {"dist":mtx2bins(dist_ref,     2.0,  20.0, 37, mask=mask),
               "omega":mtx2bins(omega_ref, -np.pi, np.pi, 25, mask=mask),
               "theta": mtx2bins(theta_ref, -np.pi, np.pi, 25, mask=(p_dist[...,0]==1)),
               "phi":mtx2bins(phi_ref,      0.0, np.pi, 13, mask=(p_dist[...,0]==1))}
      return feats

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
