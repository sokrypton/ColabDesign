import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.shared.utils import copy_dict, update_dict
from colabdesign.shared.prep import rewire
from colabdesign.af.alphafold.common import residue_constants

aa_order = residue_constants.restype_order
order_aa = {b:a for a,b in aa_order.items()}

class design_model:
  def set_opt(self, *args, **kwargs):
    '''
    set [opt]ions
    -------------------
    note: model.restart() resets the [opt]ions to their defaults
    use model.set_opt(..., set_defaults=True) 
    or model.restart(..., reset_opt=False) to avoid this
    -------------------    
    model.set_opt(models=1, recycles=0)
    model.set_opt(con=dict(num=1)) or set_opt({"con":{"num":1}})
    model.set_opt(lr=1, set_defaults=True)
    '''
    if "optimizer" in kwargs:
      print("ERROR: use model.restart(optimizer=...) to set the optimizer")

    if "recycle_mode" in kwargs and "recycle_mode" in self._args:
      if kwargs["recycle_mode"] in ["sample","last"] and self._args["recycle_mode"] in ["sample","last"]:
        self._args["recycle_mode"] = kwargs.pop("recycle_mode")
      else:
        print(f"ERROR: use {self.__class__.__name__}(recycle_mode=...) to set the recycle_mode")

    if kwargs.pop("set_defaults", False):
      update_dict(self._opt, *args, **kwargs)

    update_dict(self.opt, *args, **kwargs)

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

  def set_seq(self, seq=None, mode=None, bias=None, rm_aa=None, **kwargs):
    '''set sequence params and bias'''

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
      b, set_bias = jnp.zeros(20), False
    else:
      b, set_bias = jnp.asarray(bias), True
    
    # disable certain amino acids
    if rm_aa is not None:
      for aa in rm_aa.split(","):
        b = b.at[...,aa_order[aa]].add(-1e7)
      set_bias = True

    # use wildtype sequence
    if "wildtype" in mode or "wt" in mode:
      wt_seq = jax.nn.one_hot(self._wt_aatype,20)
      if "pos" in self.opt:
        seq = jnp.zeros(shape).at[...,self.opt["pos"],:].set(wt_seq)
      else:
        seq = wt_seq
    
    # initlaize sequence
    if seq is None:
      if hasattr(self,"key"):
        x = 0.01 * jax.random.normal(self.key(),shape)
      else:
        x = jnp.zeros(shape)
    else:
      if isinstance(seq, str):
        seq = [seq]
      if isinstance(seq, list):
        if isinstance(seq[0], str):
          seq = jnp.asarray([[aa_order.get(aa,-1) for aa in s] for s in seq])
          seq = jax.nn.one_hot(seq,20)
        else:
          seq = jnp.asarray(seq)
      
      if kwargs.pop("add_seq",False):
        b = b + seq * 1e7
        set_bias = True
      
      x = jnp.broadcast_to(jnp.asarray(seq),shape)

    if "gumbel" in mode:
      x_gumbel = jax.random.gumbel(self.key(),shape)
      if "soft" in mode:
        x = jax.nn.softmax(x + b + x_gumbel)
      elif "alpha" in self.opt:
        x = x + x_gumbel / self.opt["alpha"]
      else:
        x = x + x_gumbel

    self.params["seq"] = x
    if set_bias: self.opt["bias"] = b 

  def get_seq(self, get_best=True):
    '''
    get sequences as strings
    - set get_best=False, to get the last sampled sequence
    '''
    aux = self.aux if (self._best_aux is None or not get_best) else self._best_aux
    aux = jax.tree_map(lambda x:np.asarray(x), aux)
    x = aux["seq"]["hard"].argmax(-1)
    return ["".join([order_aa[a] for a in s]) for s in x]
  
  def get_seqs(self, get_best=True):
    return self.get_seq(get_best)

  def rewire(self, order=None, offset=0, loops=0):
    self.opt["pos"] = rewire(length=self._pos_info["length"], order=order,
                             offset=offset, loops=loops)
    if hasattr(self,"_opt"):
      self._opt["pos"] = self.opt["pos"]

  def predict(self, seq=None, model=0):
    self.set_seq(seq)
    self.set_opt(dropout=False)
    self.run(model=model, backprop=False)
    return self.aux

def soft_seq(x, opt):
  seq = {"input":x}
  # shuffle msa (randomly pick which sequence is query)
  if x.ndim == 3 and x.shape[0] > 1:
    n = jax.random.randint(key,[],0,x.shape[0])
    seq["input"] = seq["input"].at[0].set(seq["input"][n]).at[n].set(seq["input"][0])

  # straight-through/reparameterization
  seq["logits"] = seq["input"] * opt["alpha"] + opt["bias"]
  seq["soft"] = jax.nn.softmax(seq["logits"] / opt["temp"])
  seq["hard"] = jax.nn.one_hot(seq["soft"].argmax(-1), 20)
  seq["hard"] = jax.lax.stop_gradient(seq["hard"] - seq["soft"]) + seq["soft"]

  # create pseudo sequence
  seq["pseudo"] = opt["soft"] * seq["soft"] + (1-opt["soft"]) * seq["input"]
  seq["pseudo"] = opt["hard"] * seq["hard"] + (1-opt["hard"]) * seq["pseudo"]
  return seq