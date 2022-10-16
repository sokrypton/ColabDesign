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
from colabdesign.shared.model import design_model, soft_seq

# borrow some stuff from AfDesign
from colabdesign.af.prep import prep_pdb, order_aa
from colabdesign.af.alphafold.common import protein, residue_constants

from scipy.special import softmax, log_softmax

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
    self._tied_lengths = False

  def prep_inputs(self, pdb_filename=None, chain=None, homooligomer=False,
                  ignore_missing=True, fix_pos=None, inverse=False,
                  rm_aa=None, **kwargs):
    
    '''get inputs from input pdb'''
    pdb = prep_pdb(pdb_filename, chain, ignore_missing=ignore_missing)
    atom_idx = tuple(residue_constants.atom_order[k] for k in ["N","CA","C","O"])
    chain_idx = np.concatenate([[n]*l for n,l in enumerate(pdb["lengths"])])
    self._inputs = {"X":           pdb["batch"]["all_atom_positions"][:,atom_idx],
                    "mask":        pdb["batch"]["all_atom_mask"][:,1],
                    "S":           pdb["batch"]["aatype"],
                    "residue_idx": pdb["residue_index"],
                    "chain_idx":   chain_idx}
    
    self._lengths = pdb["lengths"]
    self._len = sum(self._lengths)
    self.set_seq(self._inputs["S"], rm_aa=rm_aa)
    
    if fix_pos is not None:
      p = prep_pos(fix_pos, **pdb["idx"])["pos"]
      if inverse:
        p = np.delete(np.arange(self._len),p)
      self._inputs["fix_pos"] = p
      self._inputs["bias"][p] = 1e7 * np.eye(21)[self._inputs["S"]][p,:20]
    
    if homooligomer:
      assert min(self._lengths) == max(self._lengths)
      self._tied_lengths = True

    self.pdb = pdb

  def get_af_inputs(self, af):
    '''get inputs from alphafold model'''

    self._inputs["residue_idx"] = af._inputs["residue_index"]
    self._inputs["chain_idx"]   = af._inputs["asym_id"]
    self._inputs["bias"]        = af._inputs["bias"]

    if "batch" in af._inputs:
      atom_idx = tuple(residue_constants.atom_order[k] for k in ["N","CA","C","O"])
      batch = af._inputs["batch"]
      self._inputs["X"]    = batch["all_atom_positions"][:,atom_idx]
      self._inputs["mask"] = batch["all_atom_mask"][:,1]
      self._inputs["S"]    = batch["aatype"]

    self._lengths = af._lengths
    self._len = af._len

    if "fix_pos" in af.opt:
      self._inputs["fix_pos"] = p = af.opt["fix_pos"]
      self._inputs["bias"][p] = 1e7 * np.eye(21)[self._inputs["S"]][p,:20]

    if af._args["homooligomer"]:
      assert min(self._lengths) == max(self._lengths)
      self._tied_lengths = True

  def sample(self, temperature=0.1, rescore=False, **kwargs):
    '''sample sequence'''
    O = self.sample_parallel(temperature, batch=1, rescore=rescore, **kwargs)
    return jax.tree_map(lambda x:x[0], O)
    
  def sample_parallel(self, temperature=0.1, batch=10, rescore=False, **kwargs):
    '''sample new sequence(s) in parallel'''
    I = copy_dict(self._inputs)
    I.update(kwargs)
    key = I.pop("key",self.key())
    keys = jax.random.split(key,batch)
    O = self._sample_parallel(keys, I, temperature, tied_lengths=self._tied_lengths)
    if rescore:
      O = self._rescore_parallel(keys, O["S"], O["decoding_order"], I)
    O = jax.tree_map(np.array, O)
    O["seq"] = _get_seq(O)
    O["score"] = _get_score(I,O)
    return O

  def score(self, seq=None, **kwargs):
    '''score sequence'''
    I = copy_dict(self._inputs)
    if seq is not None:
      self.set_seq(seq)
      I["S"] = self._params["seq"][0]
    I.update(kwargs)
    key = I.pop("key",self.key())
    O = jax.tree_map(np.array, self._score(**I, key=key))
    O["score"] = _get_score(I,O)
    return O

  def get_logits(self, **kwargs):
    '''get logits'''
    return self.score(**kwargs)["logits"]

  def get_unconditional_logits(self, **kwargs):
    kwargs["decoding_order"] = np.full(self._len,-1)
    return self.score(**kwargs)["logits"]

  def _setup(self):
    def _score(X, mask, residue_idx, chain_idx, key, **kwargs):
      I = {'X': X,
           'mask': mask,
           'residue_idx': residue_idx,
           'chain_idx': chain_idx}
      I.update(kwargs)
      for k in ["S","bias"]:
        if k in I: I[k] = _aa_convert(I[k])

      O = self._model.score(self._model.params, key, I)
      O["S"] = _aa_convert(O["S"], rev=True)
      O["logits"] = _aa_convert(O["logits"], rev=True)
      return O
    
    def _sample(X, mask, residue_idx, chain_idx, key, temperature=0.1, **kwargs):
      
      I = {'X': X,
           'mask': mask,
           'residue_idx': residue_idx,
           'chain_idx': chain_idx,
           'temperature': temperature}

      I.update(kwargs)
      for k in ["S","bias"]:
        if k in I: I[k] = _aa_convert(I[k])
      
      O = self._model.sample(self._model.params, key, I)
      O["S"] = _aa_convert(O["S"], rev=True)
      O["logits"] = _aa_convert(O["logits"], rev=True)
      return O

    self._score = jax.jit(_score)
    self._sample = jax.jit(_sample, static_argnames="tied_lengths")

    def _sample_parallel(key, inputs, temp, tied_lengths=False):
      inputs.pop("temperature",None)
      inputs.pop("key",None)
      return _sample(**inputs, temperature=temp, key=key, tied_lengths=tied_lengths)
    fn = jax.vmap(_sample_parallel, in_axes=[0,None,None,None])
    self._sample_parallel = jax.jit(fn, static_argnames="tied_lengths")

    def _rescore_parallel(key, S, decoding_order, inputs):
      inputs.pop("S",None)
      inputs.pop("decoding_order",None)
      inputs.pop("key",None)
      return _score(**inputs, S=S, decoding_order=decoding_order, key=key)
    fn = jax.vmap(_rescore_parallel, in_axes=[0,0,0,None])
    self._rescore_parallel = jax.jit(fn)

#######################################################################################

def _get_seq(O):
  seqs, S = [], O["S"]
  for seq in S.argmax(-1):
    seqs.append("".join([order_aa[a] for a in seq]))
  return np.array(seqs)

def _get_score(I,O):
  log_q = log_softmax(O["logits"],-1)[...,:20]
  q = softmax(O["logits"][...,:20],-1)
  if "S" in O:
    p = O["S"][...,:20]
    score = -(p * log_q).sum(-1)
  else:
    score = -(q * log_q).sum(-1)
    
  mask = I["mask"].copy()
  if "fix_pos" in I:
    mask[I["fix_pos"]] = 0
  score = (score * mask).sum(-1) / mask.sum()
  return score

def _aa_convert(x, rev=False):
  mpnn_alphabet = 'ACDEFGHIKLMNPQRSTVWYX'
  af_alphabet =   'ARNDCQEGHILKMFPSTWYVX'
  if x is None:
    return x
  else:
    if rev:
      return x[...,tuple(mpnn_alphabet.index(k) for k in af_alphabet)]
    else:
      x = jax.nn.one_hot(x,21) if jnp.issubdtype(x.dtype, jnp.integer) else x
      if x.shape[-1] == 20:
        x = jnp.pad(x,[[0,0],[0,1]])
      return x[...,tuple(af_alphabet.index(k) for k in mpnn_alphabet)]