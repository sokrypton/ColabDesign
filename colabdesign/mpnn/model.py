import jax
import jax.numpy as jnp
import numpy as np
import re
import copy
import random
import os
import joblib

from .weights import __file__ as mpnn_path
from .modules import RunModel

from colabdesign.shared.prep import prep_pos
from colabdesign.shared.utils import Key, copy_dict

# borrow some stuff from AfDesign
from colabdesign.af.prep import prep_pdb
from colabdesign.af.alphafold.common import protein, residue_constants
aa_order = residue_constants.restype_order
order_aa = {b:a for a,b in aa_order.items()}

from scipy.special import softmax, log_softmax

class mk_mpnn_model():
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
    self._inputs = {}
    self._tied_lengths = False

  def prep_inputs(self, pdb_filename=None, chain=None, homooligomer=False,
                  ignore_missing=True, fix_pos=None, inverse=False,
                  rm_aa=None, verbose=False, **kwargs):
    
    '''get inputs from input pdb'''
    pdb = prep_pdb(pdb_filename, chain, ignore_missing=ignore_missing)
    atom_idx = tuple(residue_constants.atom_order[k] for k in ["N","CA","C","O"])
    chain_idx = np.concatenate([[n]*l for n,l in enumerate(pdb["lengths"])])
    self._lengths = pdb["lengths"]
    L = sum(self._lengths)

    self._inputs = {"X":           pdb["batch"]["all_atom_positions"][:,atom_idx],
                    "mask":        pdb["batch"]["all_atom_mask"][:,1],
                    "S":           pdb["batch"]["aatype"],
                    "residue_idx": pdb["residue_index"],
                    "chain_idx":   chain_idx,
                    "lengths":     np.array(self._lengths),
                    "bias":        np.zeros((L,20))}
    

    if rm_aa is not None:
      for aa in rm_aa.split(","):
        self._inputs["bias"][...,aa_order[aa]] -= 1e6
    
    if fix_pos is not None:
      p = prep_pos(fix_pos, **pdb["idx"])["pos"]
      if inverse:
        p = np.delete(np.arange(L),p)
      self._inputs["fix_pos"] = p
      self._inputs["bias"][p] = 1e7 * np.eye(21)[self._inputs["S"]][p,:20]
    
    if homooligomer:
      assert min(self._lengths) == max(self._lengths)
      self._tied_lengths = True
      self._len = self._lengths[0]
    else:
      self._tied_lengths = False    
      self._len = sum(self._lengths)  

    self.pdb = pdb
    
    if verbose:
      print("lengths", self._lengths)
      if "fix_pos" in self._inputs:
        print("the following positions will be fixed:")
        print(self._inputs["fix_pos"])

  def get_af_inputs(self, af):
    '''get inputs from alphafold model'''

    self._lengths = af._lengths
    self._len = af._len

    self._inputs["residue_idx"] = af._inputs["residue_index"]
    self._inputs["chain_idx"]   = af._inputs["asym_id"]
    self._inputs["lengths"]     = np.array(self._lengths)

    # set bias
    L = sum(self._lengths)
    self._inputs["bias"] = np.zeros((L,20))
    self._inputs["bias"][-af._len:] = af._inputs["bias"]
    
    if "offset" in af._inputs:
      self._inputs["offset"] = af._inputs["offset"]

    if "batch" in af._inputs:
      atom_idx = tuple(residue_constants.atom_order[k] for k in ["N","CA","C","O"])
      batch = af._inputs["batch"]
      self._inputs["X"]    = batch["all_atom_positions"][:,atom_idx]
      self._inputs["mask"] = batch["all_atom_mask"][:,1]
      self._inputs["S"]    = batch["aatype"]

    # fix positions
    if af.protocol == "binder":
      p = np.arange(af._target_len)
    else:
      p = af.opt.get("fix_pos",None)
    
    if p is not None:
      self._inputs["fix_pos"] = p
      self._inputs["bias"][p] = 1e7 * np.eye(21)[self._inputs["S"]][p,:20]

    # tie positions
    if af._args["homooligomer"]:
      assert min(self._lengths) == max(self._lengths)
      self._tied_lengths = True
    else:
      self._tied_lengths = False

  def sample(self, num=1, batch=1, temperature=0.1, rescore=False, **kwargs):
    '''sample sequence'''
    O = []
    for _ in range(num):
      O.append(self.sample_parallel(batch, temperature, rescore, **kwargs))
    return jax.tree_map(lambda *x:np.concatenate(x,0),*O)    

  def sample_parallel(self, batch=10, temperature=0.1, rescore=False, **kwargs):
    '''sample new sequence(s) in parallel'''
    I = copy_dict(self._inputs)
    I.update(kwargs)
    key = I.pop("key",self.key())
    keys = jax.random.split(key,batch)
    O = self._sample_parallel(keys, I, temperature, self._tied_lengths)
    if rescore:
      O = self._rescore_parallel(keys, I, O["S"], O["decoding_order"])
    O = jax.tree_map(np.array, O)

    # process outputs to human-readable form
    O.update(self._get_seq(O))
    O.update(self._get_score(I,O))
    return O

  def _get_seq(self, O):
    ''' one_hot to amino acid sequence '''
    def split_seq(seq):
      if len(self._lengths) > 1:
        seq = "".join(np.insert(list(seq),np.cumsum(self._lengths[:-1]),"/"))
        if self._tied_lengths:
          seq = seq.split("/")[0]
      return seq
    seqs, S = [], O["S"].argmax(-1)
    if S.ndim == 1: S = [S]
    for s in S:
      seq = "".join([order_aa[a] for a in s])
      seq = split_seq(seq)
      seqs.append(seq)
    
    return {"seq": np.array(seqs)}

  def _get_score(self, I, O):
    ''' logits to score/sequence_recovery '''
    mask = I["mask"].copy()
    if "fix_pos" in I:
      mask[I["fix_pos"]] = 0

    log_q = log_softmax(O["logits"],-1)[...,:20]
    q = softmax(O["logits"][...,:20],-1)
    if "S" in O:
      S = O["S"][...,:20]
      score = -(S * log_q).sum(-1)
      seqid = S.argmax(-1) == self._inputs["S"]
    else:
      score = -(q * log_q).sum(-1)
      seqid = np.zeros_like(score)
      
    score = (score * mask).sum(-1) / (mask.sum() + 1e-8)
    seqid = (seqid * mask).sum(-1) / (mask.sum() + 1e-8)

    return {"score":score, "seqid":seqid}

  def score(self, seq=None, **kwargs):
    '''score sequence'''
    I = copy_dict(self._inputs)
    if seq is not None:
      p = np.arange(I["S"].shape[0])
      if self._tied_lengths and len(seq) == self._lengths[0]:
        seq = seq * len(self._lengths)
      if "fix_pos" in I and len(seq) == (I["S"].shape[0] - I["fix_pos"].shape[0]):
        p = np.delete(p,I["fix_pos"])
      I["S"][p] = np.array([aa_order.get(aa,-1) for aa in seq])
    I.update(kwargs)
    key = I.pop("key",self.key())
    O = jax.tree_map(np.array, self._score(**I, key=key))
    O.update(self._get_score(I,O))
    return O

  def get_logits(self, **kwargs):
    '''get logits'''
    return self.score(**kwargs)["logits"]

  def get_unconditional_logits(self, **kwargs):
    L = self._inputs["X"].shape[0]
    kwargs["ar_mask"] = np.zeros((L,L))
    return self.score(**kwargs)["logits"]

  def set_seed(self, seed=None):
    np.random.seed(seed=seed)
    self.key = Key(seed=seed).get

  def _setup(self):
    def _score(X, mask, residue_idx, chain_idx, key, **kwargs):
      I = {'X': X,
           'mask': mask,
           'residue_idx': residue_idx,
           'chain_idx': chain_idx}
      I.update(kwargs)

      # define decoding order
      if "decoding_order" not in I:
        key, sub_key = jax.random.split(key)
        randn = jax.random.uniform(sub_key, (I["X"].shape[0],))    
        randn = jnp.where(I["mask"], randn, randn+1)
        if "fix_pos" in I: randn = randn.at[I["fix_pos"]].add(-1)        
        I["decoding_order"] = randn.argsort()

      for k in ["S","bias"]:
        if k in I: I[k] = _aa_convert(I[k])

      O = self._model.score(self._model.params, key, I)
      O["S"] = _aa_convert(O["S"], rev=True)
      O["logits"] = _aa_convert(O["logits"], rev=True)
      return O
    
    def _sample(X, mask, residue_idx, chain_idx, key,
                temperature=0.1, tied_lengths=False, **kwargs):
      I = {'X': X,
           'mask': mask,
           'residue_idx': residue_idx,
           'chain_idx': chain_idx,
           'temperature': temperature}
      I.update(kwargs)

      # define decoding order
      if "decoding_order" in I:
        if I["decoding_order"].ndim == 1:
          I["decoding_order"] = I["decoding_order"][:,None]
      else:
        key, sub_key = jax.random.split(key)
        randn = jax.random.uniform(sub_key, (I["X"].shape[0],))    
        randn = jnp.where(I["mask"], randn, randn+1)
        if "fix_pos" in I: randn = randn.at[I["fix_pos"]].add(-1)        
        if tied_lengths:
          copies = I["lengths"].shape[0]
          decoding_order_tied = randn.reshape(copies,-1).mean(0).argsort()
          I["decoding_order"] = jnp.arange(I["X"].shape[0]).reshape(copies,-1).T[decoding_order_tied]
        else:
          I["decoding_order"] = randn.argsort()[:,None]

      for k in ["S","bias"]:
        if k in I: I[k] = _aa_convert(I[k])
      
      O = self._model.sample(self._model.params, key, I)
      O["S"] = _aa_convert(O["S"], rev=True)
      O["logits"] = _aa_convert(O["logits"], rev=True)
      return O

    self._score = jax.jit(_score)
    self._sample = jax.jit(_sample, static_argnames=["tied_lengths"])

    def _sample_parallel(key, inputs, temperature, tied_lengths=False):
      inputs.pop("temperature",None)
      inputs.pop("key",None)
      return _sample(**inputs, key=key, temperature=temperature, tied_lengths=tied_lengths)
    fn = jax.vmap(_sample_parallel, in_axes=[0,None,None,None])
    self._sample_parallel = jax.jit(fn, static_argnames=["tied_lengths"])

    def _rescore_parallel(key, inputs, S, decoding_order):
      inputs.pop("S",None)
      inputs.pop("decoding_order",None)
      inputs.pop("key",None)
      return _score(**inputs, key=key, S=S, decoding_order=decoding_order)
    fn = jax.vmap(_rescore_parallel, in_axes=[0,None,0,0])
    self._rescore_parallel = jax.jit(fn)

#######################################################################################

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