import jax
import jax.numpy as jnp
import numpy as np
import re
import copy
import random
import os
import joblib

from .jax_weights import __file__ as mpnn_path
from .modules import RunModel

from colabdesign.shared.prep import prep_pos
from colabdesign.shared.utils import Key, copy_dict

# borrow some stuff from AfDesign
from colabdesign.af.prep import prep_pdb, order_aa
from colabdesign.af.alphafold.common import protein, residue_constants
from colabdesign.shared.model import design_model, soft_seq

class mk_mpnn_model(design_model):
  def __init__(self, model_name="v_48_020",
               backbone_noise=0.0, dropout=0.0,
               seed=None, verbose=False):
    # load model
    path = os.path.join(os.path.dirname(mpnn_path), f'{model_name}.pkl')    
    checkpoint = joblib.load(path)
    config = {'num_letters': 21,
              'node_features': 128,
              'edge_features': 128,
              'hidden_dim': 128,
              'num_encoder_layers': 3,
              'num_decoder_layers': 3,
              'augment_eps': backbone_noise,
              'k_neighbors': checkpoint['num_edges'],
              'dropout': dropout}
    self._model = RunModel(config)
    self._model.params = jax.tree_map(np.array, checkpoint['model_state_dict'])
    self._setup()
    self.set_seed(seed)

    self._num = 1
    self._params = {}
    self._inputs = {}

  def prep_inputs(self, pdb_filename=None, chain=None, fix_pos=None, ignore_missing=True,
                  rm_aa=None, **kwargs):
    '''get inputs from input pdb'''
    pdb = prep_pdb(pdb_filename, chain, ignore_missing=ignore_missing)
    atom_idx = tuple(residue_constants.atom_order[k] for k in ["N","CA","C","O"])
    chain_idx = np.concatenate([[n]*l for n,l in enumerate(pdb["lengths"])])
    self._inputs = {"X":           pdb["batch"]["all_atom_positions"][:,atom_idx],
                    "S":           pdb["batch"]["aatype"],
                    "mask":        pdb["batch"]["all_atom_mask"][:,1],
                    "residue_idx": pdb["residue_index"],
                    "chain_idx":   chain_idx}
    self._len = self._inputs["X"].shape[0]
    self.set_seq(self._inputs["S"], rm_aa=rm_aa)
    if fix_pos is not None:
      self._inputs["fix_pos"] = p = prep_pos(fix_pos, **pdb["idx"])["pos"]
      self._inputs["bias"][p] = 1e7 * np.eye(21)[self._inputs["S"]][p,:20]

  def get_af_inputs(self, af):
    '''get inputs from alphafold model'''
    atom_idx = tuple(residue_constants.atom_order[k] for k in ["N","CA","C","O"])
    batch = af._inputs["batch"]
    self._inputs = {"X":           batch["all_atom_positions"][:,atom_idx],
                    "S":           batch["aatype"],
                    "mask":        batch["all_atom_mask"][:,1],
                    "residue_idx": af._inputs["residue_index"],
                    "chain_idx":   af._inputs["asym_id"],
                    "bias":        af._inputs["bias"]}
    self._len = self._inputs["X"].shape[0]
    if "fix_pos" in af.opt:
      self._inputs["fix_pos"] = p = af.opt["fix_pos"]
      self._inputs["bias"][p] = 1e7 * np.eye(21)[self._inputs["S"]][p,:20]
  
  def sample(self, temperature=0.1, **kwargs):
    '''sample sequence'''
    self.sample_parallel(temperature, batch=1, **kwargs)
    self._outputs = jax.tree_map(lambda x:x[0], self._outputs)
    return self._outputs
    
  def sample_parallel(self, temperature=0.1, batch=10, **kwargs):
    '''sample new sequence(s) in parallel'''

    # define parallel function
    if not hasattr(self,"_sample_parallel"):
      def _sample(key, inputs, temp):
        return self._sample(**inputs, temperature=temp, key=key)
      self._sample_parallel = jax.jit(jax.vmap(_sample, in_axes=[0,None,None]))

    inputs = copy_dict(self._inputs)
    inputs.update(kwargs)
    if "fix_pos" in inputs and "decoding_order" not in inputs:
      i, p = np.arange(inputs["X"].shape[0]), inputs["fix_pos"]
      inputs["decoding_order"] = np.append(np.take(i,p),np.random.permutation(np.delete(i,p)))
    key = inputs.pop("key",self.key())
    O = self._sample_parallel(jax.random.split(key,batch), inputs, temperature)
    O = jax.tree_map(np.array, O)
    O["seq"] = np.array(["".join([order_aa[a] for a in x]) for x in O["S"].argmax(-1)])
    self._outputs = O
    return O

  def score(self, seq=None, **kwargs):
    '''score sequence'''
    inputs = copy_dict(self._inputs)
    if seq is not None:
      self.set_seq(seq)
      inputs["S"] = self._params["seq"][0]
    inputs.update(kwargs)
    if "fix_pos" in inputs and "decoding_order" not in inputs:
      i, p = np.arange(inputs["X"].shape[0]), inputs["fix_pos"]
      inputs["decoding_order"] = np.append(np.take(i,p),np.random.permutation(np.delete(i,p)))
    key = inputs.pop("key",self.key())
    return jax.tree_map(np.array, self._score(**inputs, key=key))

  def get_logits(self, **kwargs):
    '''get logits'''
    return self.score(**kwargs)["logits"]

  def get_unconditional_logits(self, **kwargs):
    '''get unconditional logits'''
    kwargs["S"], kwargs["seq"] = None, None
    return self.score(**kwargs)["logits"]

  def _setup(self):
    def _aa_convert(x, rev=False):
      mpnn_alphabet = 'ACDEFGHIKLMNPQRSTVWYX'
      af_alphabet =   'ARNDCQEGHILKMFPSTWYVX'
      if rev:
        return x[...,tuple(mpnn_alphabet.index(k) for k in af_alphabet)]
      else:
        if x.shape[-1] == 20:
          x = jnp.pad(x,[[0,0],[0,1]])
        return x[...,tuple(af_alphabet.index(k) for k in mpnn_alphabet)]

    def _score(X, mask, residue_idx, chain_idx, key, S=None,
               ar_mask=None, decoding_order=None, offset=None,
               fix_pos=None, **kwargs):
      I = {'X': X, 
           'mask': mask,
           'residue_idx': residue_idx,
           'chain_idx': chain_idx,
           'ar_mask': ar_mask,
           'offset': offset,
           'decoding_order':decoding_order}
      if S is not None:
        one_hot = jax.nn.one_hot(S,21) if jnp.issubdtype(S.dtype, jnp.integer) else S
        I["S"] = _aa_convert(one_hot)

      O = self._model.score(self._model.params, key, I)
      O["logits"] = _aa_convert(O["logits"], rev=True)
      if S is not None:
        seq = _aa_convert(I["S"], rev=True)
        O["score"] = -(seq * jax.nn.log_softmax(O["logits"])).sum(-1)
      else:
        O["score"] = -jax.nn.log_softmax(O["logits"])[:,:20].max(-1)

      if fix_pos is not None: mask = mask.at[fix_pos].set(0)
      O["score"] = (O["score"] * mask).sum()/mask.sum()
      return O
    
    def _sample(X, mask, residue_idx, chain_idx, key,
                decoding_order=None, offset=None, temperature=0.1,
                bias=None, fix_pos=None, **kwargs):
      # sample input
      I = {'X': X,
           'mask': mask,
           'residue_idx': residue_idx,
           'chain_idx': chain_idx,
           'temperature': temperature,
           'decoding_order': decoding_order,
           'offset': offset}
      if bias is not None:
        I["bias"] = _aa_convert(bias)

      O = self._model.sample(self._model.params, key, I)
      O["S"] = _aa_convert(O["S"], rev=True)
      O["logits"] = _aa_convert(O["logits"], rev=True)
      O["score"] = -(O["S"] * jax.nn.log_softmax(O["logits"])).sum(-1)

      if fix_pos is not None: mask = mask.at[fix_pos].set(0)
      O["score"] = (O["score"] * mask).sum()/mask.sum()
      return O

    ###########################################################################

    self._score = jax.jit(_score)
    self._sample = jax.jit(_sample)