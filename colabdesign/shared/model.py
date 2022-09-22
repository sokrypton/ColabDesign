import jax
import jax.numpy as jnp
import numpy as np
import optax

from colabdesign.shared.utils import copy_dict, update_dict, softmax, Key
from colabdesign.shared.prep import rewire
from colabdesign.af.alphafold.common import residue_constants

aa_order = residue_constants.restype_order
order_aa = {b:a for a,b in aa_order.items()}

class design_model:
  def set_weights(self, *args, **kwargs):
    '''
    set weights
    -------------------
    note: model.restart() resets the weights to their defaults
    use model.set_weights(..., set_defaults=True) to avoid this
    -------------------
    model.set_weights(rmsd=1)
    '''
    if kwargs.pop("set_defaults", False):
      update_dict(self._opt["weights"], *args, **kwargs)
    update_dict(self.opt["weights"], *args, **kwargs)

  def set_seq(self, seq=None, mode=None, bias=None, rm_aa=None, set_state=True, **kwargs):
    '''
    set sequence params and bias
    -----------------------------------
    -seq=str or seq=[str,str] or seq=array(shape=(L,20) or shape=(?,L,20))
    -mode=
      -"wildtype"/"wt" = initialize sequence with sequence saved from input PDB
      -"gumbel" = initial sequence with gumbel distribution
      -"soft_???" = apply softmax-activation to initailized sequence (eg. "soft_gumbel")
    -bias=array(shape=(20,) or shape=(L,20)) - bias the sequence
    -rm_aa="C,W" = specify which amino acids to remove (aka. add a negative-infinity bias to these aa)
    -----------------------------------
    '''

    # backward compatibility
    seq_init = kwargs.pop("seq_init",None)
    if seq_init is not None:
      modes = ["soft","gumbel","wildtype","wt"]
      if isinstance(seq_init,str):
        seq_init = seq_init.split("_")
      if isinstance(seq_init,list) and seq_init[0] in modes:
        mode = seq_init
      else:
        seq = seq_init

    if mode is None: mode = []

    # decide on shape
    shape = (self._num, self._len, 20)
    
    # initialize bias
    if bias is None:
      b, set_bias = np.zeros(20), False
    else:
      b, set_bias = np.asarray(bias), True
    
    # disable certain amino acids
    if rm_aa is not None:
      for aa in rm_aa.split(","):
        b = b.at[...,aa_order[aa]].add(-1e6)
      set_bias = True

    # use wildtype sequence
    if ("wildtype" in mode or "wt" in mode) and hasattr(self,"_wt_aatype"):
      wt_seq = np.eye(20)[self._wt_aatype]
      if "pos" in self.opt and self.opt["pos"].shape[0] == wt_seq.shape[0]:
        seq = np.zeros(shape).at[...,self.opt["pos"],:].set(wt_seq)
      else:
        seq = wt_seq
    
    # initialize sequence
    if seq is None:
      if hasattr(self,"key"):
        x = 0.01 * np.asarray(jax.random.normal(self.key(),shape))
      else:
        x = np.zeros(shape)
    else:
      if isinstance(seq, str):
        seq = [seq]
      if isinstance(seq, list):
        if isinstance(seq[0], str):
          seq = np.asarray([[aa_order.get(aa,-1) for aa in s] for s in seq])
        else:
          seq = np.asarray(seq)
      else:
        seq = np.asarray(seq)

      if np.issubdtype(seq.dtype, np.integer):
        seq = np.eye(20)[seq]
      
      if kwargs.pop("add_seq",False):
        b = b + seq * 1e7
        set_bias = True
      
      x = np.broadcast_to(seq,shape)

    if "gumbel" in mode:
      y_gumbel = np.asarray(jax.random.gumbel(self.key(),shape))
      if "soft" in mode:
        y = softmax(x + b + y_gumbel)
      elif "alpha" in self.opt:
        y = x + y_gumbel / self.opt["alpha"]
      else:
        y = x + y_gumbel
      
      x = np.where(x.sum(-1,keepdims=True) == 1, x, y)      

    # set seq/bias/state
    self._params["seq"] = x
    
    if set_bias:
      self.opt["bias"] = b 

  def _norm_seq_grad(self):
    g = self.aux["grad"]["seq"]
    eff_L = (np.square(g).sum(-1,keepdims=True) > 0).sum(-2,keepdims=True)
    gn = np.linalg.norm(g,axis=(-1,-2),keepdims=True)
    self.aux["grad"]["seq"] = g * np.sqrt(eff_L) / (gn + 1e-7)  

  def set_optimizer(self, optimizer=None, learning_rate=None, **kwargs):
    '''
    set/reset optimizer
    ----------------------------------
    supported optimizers include: [adabelief, adafactor, adagrad, adam, adamw, 
    fromage, lamb, lars, noisy_sgd, dpsgd, radam, rmsprop, sgd, sm3, yogi]
    '''
    optimizers = {'adabelief':optax.adabelief,'adafactor':optax.adafactor,
                  'adagrad':optax.adagrad,'adam':optax.adam,
                  'adamw':optax.adamw,'fromage':optax.fromage,
                  'lamb':optax.lamb,'lars':optax.lars,
                  'noisy_sgd':optax.noisy_sgd,'dpsgd':optax.dpsgd,
                  'radam':optax.radam,'rmsprop':optax.rmsprop,
                  'sgd':optax.sgd,'sm3':optax.sm3,'yogi':optax.yogi}
    
    if optimizer is not None: self._args["optimizer"] = optimizer
    o = optimizers[self._args["optimizer"]](1.0, **kwargs)
    self._state = o.init(self._params)
    self.opt["learning_rate"] = 0.1 if learning_rate is None else learning_rate

    def update_grad(state, grad):
      updates, state = o.update(grad, state)
      grad = jax.tree_map(lambda x:-x, updates)
      return state, grad
    self._optimizer = jax.jit(update_grad, backend='cpu')

  def set_seed(self, seed=None):
    self.key = Key(seed=seed).get
    
  def get_seq(self, get_best=True):
    '''
    get sequences as strings
    - set get_best=False, to get the last sampled sequence
    '''
    aux = self._best["aux"] if (get_best and "aux" in self._best) else self.aux
    aux = jax.tree_map(lambda x:np.asarray(x), aux)
    x = aux["seq"]["hard"].argmax(-1)
    return ["".join([order_aa[a] for a in s]) for s in x]
  
  def get_seqs(self, get_best=True):
    return self.get_seq(get_best)

  def rewire(self, order=None, offset=0, loops=0):
    '''
    helper function for "partial" protocol
    -----------------------------------------
    -order=[0,1,2] - change order of specified segments
    -offset=0 - specify start position of the first segment
    -loops=[3,2] - specified loop lengths between segments
    -----------------------------------------
    '''
    self.opt["pos"] = rewire(length=self._pos_info["length"], order=order,
                             offset=offset, loops=loops)

    # make default
    if hasattr(self,"_opt"): self._opt["pos"] = self.opt["pos"]

def soft_seq(x, opt, key=None):
  seq = {"input":x}
  # shuffle msa (randomly pick which sequence is query)
  if x.ndim == 3 and x.shape[0] > 1 and key is not None:
    key, sub_key = jax.random.split(key)
    n = jax.random.randint(sub_key,[],0,x.shape[0])
    seq["input"] = seq["input"].at[0].set(seq["input"][n]).at[n].set(seq["input"][0])

  # straight-through/reparameterization
  seq["logits"] = seq["input"] * opt["alpha"] + opt["bias"]
  seq["pssm"] = jax.nn.softmax(seq["logits"])
  seq["soft"] = jax.nn.softmax(seq["logits"] / opt["temp"])
  seq["hard"] = jax.nn.one_hot(seq["soft"].argmax(-1), 20)
  seq["hard"] = jax.lax.stop_gradient(seq["hard"] - seq["soft"]) + seq["soft"]

  # create pseudo sequence
  seq["pseudo"] = opt["soft"] * seq["soft"] + (1-opt["soft"]) * seq["input"]
  seq["pseudo"] = opt["hard"] * seq["hard"] + (1-opt["hard"]) * seq["pseudo"]
  return seq