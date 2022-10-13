import jax
import jax.numpy as jnp
import numpy as np
import re
import copy
import random
import os
import joblib

from colabdesign.shared.prng import SafeKey
from colabdesign.mpnn.jax_weights import __file__ as mpnn_path
from .modules import RunModel

class mk_mpnn_model:
  def __init__(self, model_name="v_48_020",
               backbone_noise=0.0, dropout=0.0,
               seed=None, verbose=False):
    '''
    backbone_noise: Standard deviation of Gaussian noise to add to backbone atoms
    '''
    hidden_dim = 128
    num_layers = 3 

    path = os.path.join(os.path.dirname(mpnn_path), f'{model_name}.pkl')    
    checkpoint = joblib.load(path)
    params = jax.tree_util.tree_map(jnp.array, checkpoint['model_state_dict'])
    if verbose:
      print('Number of edges:', checkpoint['num_edges'])
      noise_level_print = checkpoint['noise_level']
      print(f'Training noise level: {noise_level_print}A')

    config = {'num_letters': 21,
              'node_features': hidden_dim,
              'edge_features': hidden_dim,
              'hidden_dim': hidden_dim,
              'num_encoder_layers': num_layers,
              'num_decoder_layers': num_layers,
              'augment_eps': backbone_noise,
              'k_neighbors': checkpoint['num_edges'],
              'dropout': dropout}
    self.model = RunModel(config)
    self.model.params = params

    if seed is None:
      seed = random.randint(0,2147483647)
    self.safe_key = SafeKey(jax.random.PRNGKey(seed))


  def _aa_convert(self, x, rev=False):
    mpnn_alphabet = 'ACDEFGHIKLMNPQRSTVWYX'
    af_alphabet = 'ARNDCQEGHILKMFPSTWYVX'
    if rev:
      return x[...,tuple(mpnn_alphabet.index(k) for k in af_alphabet)][...,:20]
    else:
      if x.shape[-1] == 20:
        x = jnp.concatenate([x,jnp.zeros_like(x[...,:1])],-1)
      return x[...,tuple(af_alphabet.index(k) for k in mpnn_alphabet)]

    
  def get_logits(self, X, mask, residue_idx, chain_idx, key, S=None,
                 ar_mask=None, decoding_order=None, offset=None, **kwargs):
    inputs = {'X': X, 'S':S, 'mask': mask,
              'residue_idx': residue_idx,
              'chain_idx': chain_idx,
              'ar_mask': ar_mask,
              'offset': offset}
    if S is not None:
      one_hot = jax.nn.one_hot(S,21) if jnp.issubdtype(S.dtype, jnp.integer) else S
      inputs["S"] = self._aa_convert(one_hot)
      if decoding_order is None:
        inputs["decoding_order"] = jnp.argsort(residue_idx,-1)
      else:
        inputs["decoding_order"] = decoding_order

    logits = self.model.score(self.model.params, key, inputs)[0]
    return self._aa_convert(logits, rev=True)

  def sample(self, X, mask, residue_idx, chain_idx, key,
             decoding_order=None, temperature=0.1, bias=None, **kwargs):

    # sample input
    inputs = {'X': X,
              'chain_idx': chain_idx,
              'residue_idx': residue_idx,
              'mask': mask,
              'temperature': temperature,
              'decoding_order': decoding_order}
    if bias is not None:
      inputs["bias"] = self._aa_convert(bias)

    sample_dict = self.model.sample(self.model.params, key, inputs)
    seq = self._aa_convert(sample_dict["S"], rev=True)
    logits = self._aa_convert(sample_dict["logits"], rev=True)
    score = -(seq * jax.nn.log_softmax(logits)).sum(-1).mean()
    return {"seq":seq, "logits":logits, "score":score}