# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Modules and code used in the core part of AlphaFold.

The structure generation code is in 'folding.py'.
"""
import functools
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.af.alphafold.model import all_atom
from colabdesign.af.alphafold.model import common_modules
from colabdesign.af.alphafold.model import folding
from colabdesign.af.alphafold.model import layer_stack
from colabdesign.af.alphafold.model import lddt
from colabdesign.af.alphafold.model import mapping
from colabdesign.af.alphafold.model import prng
from colabdesign.af.alphafold.model import quat_affine
from colabdesign.af.alphafold.model import utils
import haiku as hk
import jax
import jax.numpy as jnp

from colabdesign.af.alphafold.model.r3 import Rigids, Rots, Vecs

def apply_dropout(*, tensor, safe_key, rate, broadcast_dim=None):
  """Applies dropout to a tensor."""
  shape = list(tensor.shape)
  if broadcast_dim is not None:
    shape[broadcast_dim] = 1
  keep_rate = 1.0 - rate
  keep = jax.random.bernoulli(safe_key.get(), keep_rate, shape=shape)
  return keep * tensor / keep_rate

def dropout_wrapper(module,
                    input_act,
                    mask,
                    safe_key,
                    global_config,
                    use_dropout,
                    output_act=None,
                    **kwargs):
  """Applies module + dropout + residual update."""
  if output_act is None:
    output_act = input_act

  gc = global_config
  residual = module(input_act, mask, **kwargs)

  if module.config.shared_dropout:
    if module.config.orientation == 'per_row':
      broadcast_dim = 0
    else:
      broadcast_dim = 1
  else:
    broadcast_dim = None

  residual = apply_dropout(tensor=residual,
                           safe_key=safe_key,
                           rate=jnp.where(use_dropout, module.config.dropout_rate, 0),
                           broadcast_dim=broadcast_dim)

  new_act = output_act + residual

  return new_act


def create_extra_msa_feature(batch):
  """Expand extra_msa into 1hot and concat with other extra msa features.

  We do this as late as possible as the one_hot extra msa can be very large.

  Arguments:
    batch: a dictionary with the following keys:
     * 'extra_msa': [N_extra_seq, N_res] MSA that wasn't selected as a cluster
       centre. Note, that this is not one-hot encoded.
     * 'extra_has_deletion': [N_extra_seq, N_res] Whether there is a deletion to
       the left of each position in the extra MSA.
     * 'extra_deletion_value': [N_extra_seq, N_res] The number of deletions to
       the left of each position in the extra MSA.

  Returns:
    Concatenated tensor of extra MSA features.
  """
  # 23 = 20 amino acids + 'X' for unknown + gap + bert mask
  msa_1hot = jax.nn.one_hot(batch['extra_msa'], 23)
  msa_feat = [msa_1hot,
              jnp.expand_dims(batch['extra_has_deletion'], axis=-1),
              jnp.expand_dims(batch['extra_deletion_value'], axis=-1)]
  return jnp.concatenate(msa_feat, axis=-1)

class AlphaFoldIteration(hk.Module):
  """A single recycling iteration of AlphaFold architecture.
  Jumper et al. (2021) Suppl. Alg. 2 "Inference" lines 3-22
  """
  def __init__(self, config, global_config, name='alphafold_iteration'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, batch, **kwargs):

    # Compute representations for each batch element and average.
    evoformer_module = EmbeddingsAndEvoformer(self.config.embeddings_and_evoformer, self.global_config)
    representations = evoformer_module(batch)    
    
    head_factory = {
          'masked_msa': MaskedMsaHead,
          'distogram': DistogramHead,
          'structure_module': folding.StructureModule,
          'predicted_lddt': PredictedLDDTHead,
          'predicted_aligned_error': PredictedAlignedErrorHead,
          'experimentally_resolved': ExperimentallyResolvedHead} 

    heads = {}
    for name, head_config in sorted(self.config.heads.items()):
      if not head_config.weight: continue
      heads[name] = head_factory[name](head_config, self.global_config)
    
    ret = {'representations':representations}
    for name, head in heads.items():
      if name in ('predicted_lddt', 'predicted_aligned_error'):
        continue
      else:
        ret[name] = head(representations, batch)
        if 'representations' in ret[name]:
          representations.update(ret[name].pop('representations'))

    for name in ('predicted_lddt', 'predicted_aligned_error'):
      ret[name] = heads[name](representations, batch)
    return ret

class AlphaFold(hk.Module):
  """AlphaFold Jumper et al. (2021) Suppl. Alg. 2 "Inference"""
  def __init__(self, config, name='alphafold'):
    super().__init__(name=name)
    self.config = config
    self.global_config = config.global_config

  def __call__(self, batch, **kwargs):
    """Run the AlphaFold model."""
    impl = AlphaFoldIteration(self.config, self.global_config)

    def get_prev(ret):
      new_prev = {
          'prev_msa_first_row': ret['representations']['msa_first_row'],
          'prev_pair': ret['representations']['pair'],
          'prev_pos': ret['structure_module']['final_atom_positions']
      }
      if self.global_config.use_dgram:
        new_prev['prev_dgram'] = ret["distogram"]["logits"]
      return new_prev    

    prev = batch.pop("prev")                      
    ret = impl(batch={**batch, **prev})
    ret["prev"] = get_prev(ret)        
    return ret

class TemplatePairStack(hk.Module):
  """Pair stack for the templates.

  Jumper et al. (2021) Suppl. Alg. 16 "TemplatePairStack"
  """

  def __init__(self, config, global_config, name='template_pair_stack'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, pair_act, pair_mask, use_dropout, safe_key=None):
    """Builds TemplatePairStack module.

    Arguments:
      pair_act: Pair activations for single template, shape [N_res, N_res, c_t].
      pair_mask: Pair mask, shape [N_res, N_res].
      safe_key: Safe key object encapsulating the random number generation key.

    Returns:
      Updated pair_act, shape [N_res, N_res, c_t].
    """

    if safe_key is None:
      safe_key = prng.SafeKey(hk.next_rng_key())

    gc = self.global_config
    c = self.config

    if not c.num_block:
      return pair_act

    def block(x):
      """One block of the template pair stack."""
      pair_act, safe_key = x

      dropout_wrapper_fn = functools.partial(
          dropout_wrapper, global_config=gc, use_dropout=use_dropout)

      safe_key, *sub_keys = safe_key.split(6)
      sub_keys = iter(sub_keys)

      pair_act = dropout_wrapper_fn(
          TriangleAttention(c.triangle_attention_starting_node, gc,
                            name='triangle_attention_starting_node'),
          pair_act,
          pair_mask,
          next(sub_keys))
      pair_act = dropout_wrapper_fn(
          TriangleAttention(c.triangle_attention_ending_node, gc,
                            name='triangle_attention_ending_node'),
          pair_act,
          pair_mask,
          next(sub_keys))
      pair_act = dropout_wrapper_fn(
          TriangleMultiplication(c.triangle_multiplication_outgoing, gc,
                                 name='triangle_multiplication_outgoing'),
          pair_act,
          pair_mask,
          next(sub_keys))
      pair_act = dropout_wrapper_fn(
          TriangleMultiplication(c.triangle_multiplication_incoming, gc,
                                 name='triangle_multiplication_incoming'),
          pair_act,
          pair_mask,
          next(sub_keys))
      pair_act = dropout_wrapper_fn(
          Transition(c.pair_transition, gc, name='pair_transition'),
          pair_act,
          pair_mask,
          next(sub_keys))

      return pair_act, safe_key

    if gc.use_remat:
      block = hk.remat(block)

    res_stack = layer_stack.layer_stack(c.num_block)(block)
    pair_act, safe_key = res_stack((pair_act, safe_key))
    return pair_act


class Transition(hk.Module):
  """Transition layer.

  Jumper et al. (2021) Suppl. Alg. 9 "MSATransition"
  Jumper et al. (2021) Suppl. Alg. 15 "PairTransition"
  """

  def __init__(self, config, global_config, name='transition_block'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, act, mask):
    """Builds Transition module.

    Arguments:
      act: A tensor of queries of size [batch_size, N_res, N_channel].
      mask: A tensor denoting the mask of size [batch_size, N_res].

    Returns:
      A float32 tensor of size [batch_size, N_res, N_channel].
    """
    _, _, nc = act.shape

    num_intermediate = int(nc * self.config.num_intermediate_factor)
    mask = jnp.expand_dims(mask, axis=-1)

    act = common_modules.LayerNorm(
        axis=[-1],
        create_scale=True,
        create_offset=True,
        name='input_layer_norm')(
            act)

    transition_module = hk.Sequential([
        common_modules.Linear(
            num_intermediate,
            initializer='relu',
            name='transition1'), jax.nn.relu,
        common_modules.Linear(
            nc,
            initializer=utils.final_init(self.global_config),
            name='transition2')
    ])

    act = mapping.inference_subbatch(
        transition_module,
        self.global_config.subbatch_size,
        batched_args=[act],
        nonbatched_args=[],
        low_memory=self.global_config.subbatch_size is not None)

    return act


def glorot_uniform():
  return hk.initializers.VarianceScaling(scale=1.0,
                                         mode='fan_avg',
                                         distribution='uniform')


class Attention(hk.Module):
  """Multihead attention."""

  def __init__(self, config, global_config, output_dim, name='attention'):
    super().__init__(name=name)

    self.config = config
    self.global_config = global_config
    self.output_dim = output_dim

  def __call__(self, q_data, m_data, bias, nonbatched_bias=None):
    """Builds Attention module.

    Arguments:
      q_data: A tensor of queries, shape [batch_size, N_queries, q_channels].
      m_data: A tensor of memories from which the keys and values are
        projected, shape [batch_size, N_keys, m_channels].
      bias: A bias for the attention, shape [batch_size, N_queries, N_keys].
      nonbatched_bias: Shared bias, shape [N_queries, N_keys].

    Returns:
      A float32 tensor of shape [batch_size, N_queries, output_dim].
    """
    # Sensible default for when the config keys are missing
    key_dim = self.config.get('key_dim', int(q_data.shape[-1]))
    value_dim = self.config.get('value_dim', int(m_data.shape[-1]))
    num_head = self.config.num_head
    assert key_dim % num_head == 0
    assert value_dim % num_head == 0
    key_dim = key_dim // num_head
    value_dim = value_dim // num_head

    q_weights = hk.get_parameter(
        'query_w', shape=(q_data.shape[-1], num_head, key_dim),
        dtype=q_data.dtype,
        init=glorot_uniform())
    k_weights = hk.get_parameter(
        'key_w', shape=(m_data.shape[-1], num_head, key_dim),
        dtype=q_data.dtype,
        init=glorot_uniform())
    v_weights = hk.get_parameter(
        'value_w', shape=(m_data.shape[-1], num_head, value_dim),
        dtype=q_data.dtype,
        init=glorot_uniform())

    q = jnp.einsum('bqa,ahc->bqhc', q_data, q_weights) * key_dim**(-0.5)
    k = jnp.einsum('bka,ahc->bkhc', m_data, k_weights)
    v = jnp.einsum('bka,ahc->bkhc', m_data, v_weights)
    logits = jnp.einsum('bqhc,bkhc->bhqk', q, k) + bias
    if nonbatched_bias is not None:
      logits += jnp.expand_dims(nonbatched_bias, axis=0)

    # patch for jax > 0.3.25
    logits = jnp.clip(logits,-1e8,1e8)

    weights = jax.nn.softmax(logits)
    weighted_avg = jnp.einsum('bhqk,bkhc->bqhc', weights, v)

    if self.global_config.zero_init:
      init = hk.initializers.Constant(0.0)
    else:
      init = glorot_uniform()

    if self.config.gating:
      gating_weights = hk.get_parameter(
          'gating_w',
          shape=(q_data.shape[-1], num_head, value_dim),
          dtype=q_data.dtype,
          init=hk.initializers.Constant(0.0))
      gating_bias = hk.get_parameter(
          'gating_b',
          shape=(num_head, value_dim),
          dtype=q_data.dtype,
          init=hk.initializers.Constant(1.0))

      gate_values = jnp.einsum('bqc, chv->bqhv', q_data,
                               gating_weights) + gating_bias

      gate_values = jax.nn.sigmoid(gate_values)

      weighted_avg *= gate_values

    o_weights = hk.get_parameter(
        'output_w', shape=(num_head, value_dim, self.output_dim),
        dtype=q_data.dtype,
        init=init)
    o_bias = hk.get_parameter('output_b', shape=(self.output_dim,),
                              dtype=q_data.dtype,
                              init=hk.initializers.Constant(0.0))

    output = jnp.einsum('bqhc,hco->bqo', weighted_avg, o_weights) + o_bias

    return output


class GlobalAttention(hk.Module):
  """Global attention.

  Jumper et al. (2021) Suppl. Alg. 19 "MSAColumnGlobalAttention" lines 2-7
  """

  def __init__(self, config, global_config, output_dim, name='attention'):
    super().__init__(name=name)

    self.config = config
    self.global_config = global_config
    self.output_dim = output_dim

  def __call__(self, q_data, m_data, q_mask, bias):
    """Builds GlobalAttention module.

    Arguments:
      q_data: A tensor of queries with size [batch_size, N_queries,
        q_channels]
      m_data: A tensor of memories from which the keys and values
        projected. Size [batch_size, N_keys, m_channels]
      q_mask: A binary mask for q_data with zeros in the padded sequence
        elements and ones otherwise. Size [batch_size, N_queries, q_channels]
        (or broadcastable to this shape).
      bias: A bias for the attention.

    Returns:
      A float32 tensor of size [batch_size, N_queries, output_dim].
    """
    # Sensible default for when the config keys are missing
    key_dim = self.config.get('key_dim', int(q_data.shape[-1]))
    value_dim = self.config.get('value_dim', int(m_data.shape[-1]))
    num_head = self.config.num_head
    assert key_dim % num_head == 0
    assert value_dim % num_head == 0
    key_dim = key_dim // num_head
    value_dim = value_dim // num_head

    q_weights = hk.get_parameter(
        'query_w', shape=(q_data.shape[-1], num_head, key_dim),
        dtype=q_data.dtype,
        init=glorot_uniform())
    k_weights = hk.get_parameter(
        'key_w', shape=(m_data.shape[-1], key_dim),
        dtype=q_data.dtype,
        init=glorot_uniform())
    v_weights = hk.get_parameter(
        'value_w', shape=(m_data.shape[-1], value_dim),
        dtype=q_data.dtype,
        init=glorot_uniform())

    v = jnp.einsum('bka,ac->bkc', m_data, v_weights)

    q_avg = utils.mask_mean(q_mask, q_data, axis=1)

    q = jnp.einsum('ba,ahc->bhc', q_avg, q_weights) * key_dim**(-0.5)
    k = jnp.einsum('bka,ac->bkc', m_data, k_weights)
    bias = (1e9 * (q_mask[:, None, :, 0] - 1.))
    logits = jnp.einsum('bhc,bkc->bhk', q, k) + bias
    weights = jax.nn.softmax(logits)
    weighted_avg = jnp.einsum('bhk,bkc->bhc', weights, v)

    if self.global_config.zero_init:
      init = hk.initializers.Constant(0.0)
    else:
      init = glorot_uniform()

    o_weights = hk.get_parameter(
        'output_w', shape=(num_head, value_dim, self.output_dim),
        dtype=q_data.dtype,
        init=init)
    o_bias = hk.get_parameter('output_b', shape=(self.output_dim,),
                              dtype=q_data.dtype,
                              init=hk.initializers.Constant(0.0))

    if self.config.gating:
      gating_weights = hk.get_parameter(
          'gating_w',
          shape=(q_data.shape[-1], num_head, value_dim),
          dtype=q_data.dtype,
          init=hk.initializers.Constant(0.0))
      gating_bias = hk.get_parameter(
          'gating_b',
          shape=(num_head, value_dim),
          dtype=q_data.dtype,
          init=hk.initializers.Constant(1.0))

      gate_values = jnp.einsum('bqc, chv->bqhv', q_data, gating_weights)
      gate_values = jax.nn.sigmoid(gate_values + gating_bias)
      weighted_avg = weighted_avg[:, None] * gate_values
      output = jnp.einsum('bqhc,hco->bqo', weighted_avg, o_weights) + o_bias
    else:
      output = jnp.einsum('bhc,hco->bo', weighted_avg, o_weights) + o_bias
      output = output[:, None]
    return output


class MSARowAttentionWithPairBias(hk.Module):
  """MSA per-row attention biased by the pair representation.

  Jumper et al. (2021) Suppl. Alg. 7 "MSARowAttentionWithPairBias"
  """

  def __init__(self, config, global_config,
               name='msa_row_attention_with_pair_bias'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self,
               msa_act,
               msa_mask,
               pair_act):
    """Builds MSARowAttentionWithPairBias module.

    Arguments:
      msa_act: [N_seq, N_res, c_m] MSA representation.
      msa_mask: [N_seq, N_res] mask of non-padded regions.
      pair_act: [N_res, N_res, c_z] pair representation.

    Returns:
      Update to msa_act, shape [N_seq, N_res, c_m].
    """
    c = self.config

    assert len(msa_act.shape) == 3
    assert len(msa_mask.shape) == 2
    assert c.orientation == 'per_row'

    bias = (1e9 * (msa_mask - 1.))[:, None, None, :]
    assert len(bias.shape) == 4

    msa_act = common_modules.LayerNorm(
        axis=[-1], create_scale=True, create_offset=True, name='query_norm')(
            msa_act)

    pair_act = common_modules.LayerNorm(
        axis=[-1],
        create_scale=True,
        create_offset=True,
        name='feat_2d_norm')(
            pair_act)

    init_factor = 1. / jnp.sqrt(int(pair_act.shape[-1]))
    weights = hk.get_parameter(
        'feat_2d_weights',
        shape=(pair_act.shape[-1], c.num_head),
        dtype=msa_act.dtype,
        init=hk.initializers.RandomNormal(stddev=init_factor))
    nonbatched_bias = jnp.einsum('qkc,ch->hqk', pair_act, weights)

    attn_mod = Attention(
        c, self.global_config, msa_act.shape[-1])
    msa_act = mapping.inference_subbatch(
        attn_mod,
        self.global_config.subbatch_size,
        batched_args=[msa_act, msa_act, bias],
        nonbatched_args=[nonbatched_bias],
        low_memory=self.global_config.subbatch_size is not None)

    return msa_act


class MSAColumnAttention(hk.Module):
  """MSA per-column attention.

  Jumper et al. (2021) Suppl. Alg. 8 "MSAColumnAttention"
  """

  def __init__(self, config, global_config, name='msa_column_attention'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self,
               msa_act,
               msa_mask):
    """Builds MSAColumnAttention module.

    Arguments:
      msa_act: [N_seq, N_res, c_m] MSA representation.
      msa_mask: [N_seq, N_res] mask of non-padded regions.

    Returns:
      Update to msa_act, shape [N_seq, N_res, c_m]
    """
    c = self.config

    assert len(msa_act.shape) == 3
    assert len(msa_mask.shape) == 2
    assert c.orientation == 'per_column'

    msa_act = jnp.swapaxes(msa_act, -2, -3)
    msa_mask = jnp.swapaxes(msa_mask, -1, -2)

    bias = (1e9 * (msa_mask - 1.))[:, None, None, :]
    assert len(bias.shape) == 4

    msa_act = common_modules.LayerNorm(
        axis=[-1], create_scale=True, create_offset=True, name='query_norm')(
            msa_act)

    attn_mod = Attention(
        c, self.global_config, msa_act.shape[-1])
    msa_act = mapping.inference_subbatch(
        attn_mod,
        self.global_config.subbatch_size,
        batched_args=[msa_act, msa_act, bias],
        nonbatched_args=[],
        low_memory=self.global_config.subbatch_size is not None)

    msa_act = jnp.swapaxes(msa_act, -2, -3)

    return msa_act


class MSAColumnGlobalAttention(hk.Module):
  """MSA per-column global attention.

  Jumper et al. (2021) Suppl. Alg. 19 "MSAColumnGlobalAttention"
  """

  def __init__(self, config, global_config, name='msa_column_global_attention'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self,
               msa_act,
               msa_mask):
    """Builds MSAColumnGlobalAttention module.

    Arguments:
      msa_act: [N_seq, N_res, c_m] MSA representation.
      msa_mask: [N_seq, N_res] mask of non-padded regions.

    Returns:
      Update to msa_act, shape [N_seq, N_res, c_m].
    """
    c = self.config

    assert len(msa_act.shape) == 3
    assert len(msa_mask.shape) == 2
    assert c.orientation == 'per_column'

    msa_act = jnp.swapaxes(msa_act, -2, -3)
    msa_mask = jnp.swapaxes(msa_mask, -1, -2)

    bias = (1e9 * (msa_mask - 1.))[:, None, None, :]
    assert len(bias.shape) == 4

    msa_act = common_modules.LayerNorm(
        axis=[-1], create_scale=True, create_offset=True, name='query_norm')(
            msa_act)

    attn_mod = GlobalAttention(
        c, self.global_config, msa_act.shape[-1],
        name='attention')
    # [N_seq, N_res, 1]
    msa_mask = jnp.expand_dims(msa_mask, axis=-1)
    msa_act = mapping.inference_subbatch(
        attn_mod,
        self.global_config.subbatch_size,
        batched_args=[msa_act, msa_act, msa_mask, bias],
        nonbatched_args=[],
        low_memory=self.global_config.subbatch_size is not None)

    msa_act = jnp.swapaxes(msa_act, -2, -3)

    return msa_act


class TriangleAttention(hk.Module):
  """Triangle Attention.

  Jumper et al. (2021) Suppl. Alg. 13 "TriangleAttentionStartingNode"
  Jumper et al. (2021) Suppl. Alg. 14 "TriangleAttentionEndingNode"
  """

  def __init__(self, config, global_config, name='triangle_attention'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, pair_act, pair_mask):
    """Builds TriangleAttention module.

    Arguments:
      pair_act: [N_res, N_res, c_z] pair activations tensor
      pair_mask: [N_res, N_res] mask of non-padded regions in the tensor.

    Returns:
      Update to pair_act, shape [N_res, N_res, c_z].
    """
    c = self.config

    assert len(pair_act.shape) == 3
    assert len(pair_mask.shape) == 2
    assert c.orientation in ['per_row', 'per_column']

    if c.orientation == 'per_column':
      pair_act = jnp.swapaxes(pair_act, -2, -3)
      pair_mask = jnp.swapaxes(pair_mask, -1, -2)

    bias = (1e9 * (pair_mask - 1.))[:, None, None, :]
    assert len(bias.shape) == 4

    pair_act = common_modules.LayerNorm(
        axis=[-1], create_scale=True, create_offset=True, name='query_norm')(
            pair_act)

    init_factor = 1. / jnp.sqrt(int(pair_act.shape[-1]))
    weights = hk.get_parameter(
        'feat_2d_weights',
        shape=(pair_act.shape[-1], c.num_head),
        dtype=pair_act.dtype,
        init=hk.initializers.RandomNormal(stddev=init_factor))
    nonbatched_bias = jnp.einsum('qkc,ch->hqk', pair_act, weights)

    attn_mod = Attention(
        c, self.global_config, pair_act.shape[-1])
    pair_act = mapping.inference_subbatch(
        attn_mod,
        self.global_config.subbatch_size,
        batched_args=[pair_act, pair_act, bias],
        nonbatched_args=[nonbatched_bias],
        low_memory=self.global_config.subbatch_size is not None)

    if c.orientation == 'per_column':
      pair_act = jnp.swapaxes(pair_act, -2, -3)

    return pair_act


class MaskedMsaHead(hk.Module):
  """Head to predict MSA at the masked locations.

  The MaskedMsaHead employs a BERT-style objective to reconstruct a masked
  version of the full MSA, based on a linear projection of
  the MSA representation.
  Jumper et al. (2021) Suppl. Sec. 1.9.9 "Masked MSA prediction"
  """

  def __init__(self, config, global_config, name='masked_msa_head'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

    if global_config.multimer_mode:
      self.num_output = len(residue_constants.restypes_with_x_and_gap)
    else:
      self.num_output = config.num_output

  def __call__(self, representations, batch):
    """Builds MaskedMsaHead module.

    Arguments:
      representations: Dictionary of representations, must contain:
        * 'msa': MSA representation, shape [N_seq, N_res, c_m].
      batch: Batch, unused.

    Returns:
      Dictionary containing:
        * 'logits': logits of shape [N_seq, N_res, N_aatype] with
            (unnormalized) log probabilies of predicted aatype at position.
    """
    del batch
    logits = common_modules.Linear(
        self.num_output,
        initializer=utils.final_init(self.global_config),
        name='logits')(
            representations['msa'])
    return dict(logits=logits)

class PredictedLDDTHead(hk.Module):
  """Head to predict the per-residue LDDT to be used as a confidence measure.

  Jumper et al. (2021) Suppl. Sec. 1.9.6 "Model confidence prediction (pLDDT)"
  Jumper et al. (2021) Suppl. Alg. 29 "predictPerResidueLDDT_Ca"
  """

  def __init__(self, config, global_config, name='predicted_lddt_head'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, representations, batch):
    """Builds ExperimentallyResolvedHead module.

    Arguments:
      representations: Dictionary of representations, must contain:
        * 'structure_module': Single representation from the structure module,
             shape [N_res, c_s].
      batch: Batch, unused.

    Returns:
      Dictionary containing :
        * 'logits': logits of shape [N_res, N_bins] with
            (unnormalized) log probabilies of binned predicted lDDT.
    """
    act = representations['structure_module']

    act = common_modules.LayerNorm(
        axis=[-1],
        create_scale=True,
        create_offset=True,
        name='input_layer_norm')(
            act)

    act = common_modules.Linear(
        self.config.num_channels,
        initializer='relu',
        name='act_0')(
            act)
    act = jax.nn.relu(act)

    act = common_modules.Linear(
        self.config.num_channels,
        initializer='relu',
        name='act_1')(
            act)
    act = jax.nn.relu(act)

    logits = common_modules.Linear(
        self.config.num_bins,
        initializer=utils.final_init(self.global_config),
        name='logits')(
            act)
    # Shape (batch_size, num_res, num_bins)
    return dict(logits=logits)


class PredictedAlignedErrorHead(hk.Module):
  """Head to predict the distance errors in the backbone alignment frames.

  Can be used to compute predicted TM-Score.
  Jumper et al. (2021) Suppl. Sec. 1.9.7 "TM-score prediction"
  """

  def __init__(self, config, global_config,
               name='predicted_aligned_error_head'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, representations, batch):
    """Builds PredictedAlignedErrorHead module.

    Arguments:
      representations: Dictionary of representations, must contain:
        * 'pair': pair representation, shape [N_res, N_res, c_z].
      batch: Batch, unused.

    Returns:
      Dictionary containing:
        * logits: logits for aligned error, shape [N_res, N_res, N_bins].
        * bin_breaks: array containing bin breaks, shape [N_bins - 1].
    """

    act = representations['pair']

    # Shape (num_res, num_res, num_bins)
    logits = common_modules.Linear(
        self.config.num_bins,
        initializer=utils.final_init(self.global_config),
        name='logits')(act)
    # Shape (num_bins,)
    breaks = jnp.linspace(
        0., self.config.max_error_bin, self.config.num_bins - 1)
    return dict(logits=logits, breaks=breaks)


class ExperimentallyResolvedHead(hk.Module):
  """Predicts if an atom is experimentally resolved in a high-res structure.

  Only trained on high-resolution X-ray crystals & cryo-EM.
  Jumper et al. (2021) Suppl. Sec. 1.9.10 '"Experimentally resolved" prediction'
  """

  def __init__(self, config, global_config,
               name='experimentally_resolved_head'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, representations, batch):
    """Builds ExperimentallyResolvedHead module.

    Arguments:
      representations: Dictionary of representations, must contain:
        * 'single': Single representation, shape [N_res, c_s].
      batch: Batch, unused.

    Returns:
      Dictionary containing:
        * 'logits': logits of shape [N_res, 37],
            log probability that an atom is resolved in atom37 representation,
            can be converted to probability by applying sigmoid.
    """
    logits = common_modules.Linear(
        37,  # atom_exists.shape[-1]
        initializer=utils.final_init(self.global_config),
        name='logits')(representations['single'])
    return dict(logits=logits)


def _layer_norm(axis=-1, name='layer_norm'):
  return common_modules.LayerNorm(
      axis=axis,
      create_scale=True,
      create_offset=True,
      eps=1e-5,
      use_fast_variance=True,
      scale_init=hk.initializers.Constant(1.),
      offset_init=hk.initializers.Constant(0.),
      param_axis=axis,
      name=name)

class TriangleMultiplication(hk.Module):
  """Triangle multiplication layer ("outgoing" or "incoming").
  Jumper et al. (2021) Suppl. Alg. 11 "TriangleMultiplicationOutgoing"
  Jumper et al. (2021) Suppl. Alg. 12 "TriangleMultiplicationIncoming"
  """

  def __init__(self, config, global_config, name='triangle_multiplication'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, left_act, left_mask):
    """Builds TriangleMultiplication module.
    Arguments:
      left_act: Pair activations, shape [N_res, N_res, c_z]
      left_mask: Pair mask, shape [N_res, N_res].
    Returns:
      Outputs, same shape/type as left_act.
    """

    if self.config.fuse_projection_weights:
      return self._fused_triangle_multiplication(left_act, left_mask)
    else:
      return self._triangle_multiplication(left_act, left_mask)

  @hk.transparent
  def _triangle_multiplication(self, left_act, left_mask):
    """Implementation of TriangleMultiplication used in AF2 and AF-M<2.3."""
    c = self.config
    gc = self.global_config

    mask = left_mask[..., None]

    act = common_modules.LayerNorm(axis=[-1], create_scale=True, create_offset=True,
                       name='layer_norm_input')(left_act)
    input_act = act

    left_projection = common_modules.Linear(
        c.num_intermediate_channel,
        name='left_projection')
    left_proj_act = mask * left_projection(act)

    right_projection = common_modules.Linear(
        c.num_intermediate_channel,
        name='right_projection')
    right_proj_act = mask * right_projection(act)

    left_gate_values = jax.nn.sigmoid(common_modules.Linear(
        c.num_intermediate_channel,
        bias_init=1.,
        initializer=utils.final_init(gc),
        name='left_gate')(act))

    right_gate_values = jax.nn.sigmoid(common_modules.Linear(
        c.num_intermediate_channel,
        bias_init=1.,
        initializer=utils.final_init(gc),
        name='right_gate')(act))

    left_proj_act *= left_gate_values
    right_proj_act *= right_gate_values

    # "Outgoing" edges equation: 'ikc,jkc->ijc'
    # "Incoming" edges equation: 'kjc,kic->ijc'
    # Note on the Suppl. Alg. 11 & 12 notation:
    # For the "outgoing" edges, a = left_proj_act and b = right_proj_act
    # For the "incoming" edges, it's swapped:
    #   b = left_proj_act and a = right_proj_act
    act = jnp.einsum(c.equation, left_proj_act, right_proj_act)

    act = common_modules.LayerNorm(
        axis=[-1],
        create_scale=True,
        create_offset=True,
        name='center_layer_norm')(
            act)

    output_channel = int(input_act.shape[-1])

    act = common_modules.Linear(
        output_channel,
        initializer=utils.final_init(gc),
        name='output_projection')(act)

    gate_values = jax.nn.sigmoid(common_modules.Linear(
        output_channel,
        bias_init=1.,
        initializer=utils.final_init(gc),
        name='gating_linear')(input_act))
    act *= gate_values

    return act

  @hk.transparent
  def _fused_triangle_multiplication(self, left_act, left_mask):
    """TriangleMultiplication with fused projection weights."""
    mask = left_mask[..., None]
    c = self.config
    gc = self.global_config

    left_act = _layer_norm(axis=-1, name='left_norm_input')(left_act)

    # Both left and right projections are fused into projection.
    projection = common_modules.Linear(
        2*c.num_intermediate_channel, name='projection')
    proj_act = mask * projection(left_act)

    # Both left + right gate are fused into gate_values.
    gate_values = common_modules.Linear(
        2 * c.num_intermediate_channel,
        name='gate',
        bias_init=1.,
        initializer=utils.final_init(gc))(left_act)
    proj_act *= jax.nn.sigmoid(gate_values)

    left_proj_act = proj_act[:, :, :c.num_intermediate_channel]
    right_proj_act = proj_act[:, :, c.num_intermediate_channel:]
    act = jnp.einsum(c.equation, left_proj_act, right_proj_act)

    act = _layer_norm(axis=-1, name='center_norm')(act)

    output_channel = int(left_act.shape[-1])

    act = common_modules.Linear(
        output_channel,
        initializer=utils.final_init(gc),
        name='output_projection')(act)

    gate_values = common_modules.Linear(
        output_channel,
        bias_init=1.,
        initializer=utils.final_init(gc),
        name='gating_linear')(left_act)
    act *= jax.nn.sigmoid(gate_values)

    return act


class DistogramHead(hk.Module):
  """Head to predict a distogram.

  Jumper et al. (2021) Suppl. Sec. 1.9.8 "Distogram prediction"
  """

  def __init__(self, config, global_config, name='distogram_head'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, representations, batch):
    """Builds DistogramHead module.

    Arguments:
      representations: Dictionary of representations, must contain:
        * 'pair': pair representation, shape [N_res, N_res, c_z].
      batch: Batch, unused.

    Returns:
      Dictionary containing:
        * logits: logits for distogram, shape [N_res, N_res, N_bins].
        * bin_breaks: array containing bin breaks, shape [N_bins - 1,].
    """
    half_logits = common_modules.Linear(
        self.config.num_bins,
        initializer=utils.final_init(self.global_config),
        name='half_logits')(
            representations['pair'])

    logits = half_logits + jnp.swapaxes(half_logits, -2, -3)
    breaks = jnp.linspace(self.config.first_break, self.config.last_break,
                          self.config.num_bins - 1)

    return dict(logits=logits, bin_edges=breaks)


class OuterProductMean(hk.Module):
  """Computes mean outer product.

  Jumper et al. (2021) Suppl. Alg. 10 "OuterProductMean"
  """

  def __init__(self,
               config,
               global_config,
               num_output_channel,
               name='outer_product_mean'):
    super().__init__(name=name)
    self.global_config = global_config
    self.config = config
    self.num_output_channel = num_output_channel

  def __call__(self, act, mask):
    """Builds OuterProductMean module.

    Arguments:
      act: MSA representation, shape [N_seq, N_res, c_m].
      mask: MSA mask, shape [N_seq, N_res].

    Returns:
      Update to pair representation, shape [N_res, N_res, c_z].
    """
    gc = self.global_config
    c = self.config

    mask = mask[..., None]
    act = common_modules.LayerNorm([-1], True, True, name='layer_norm_input')(act)

    left_act = mask * common_modules.Linear(
        c.num_outer_channel,
        initializer='linear',
        name='left_projection')(
            act)

    right_act = mask * common_modules.Linear(
        c.num_outer_channel,
        initializer='linear',
        name='right_projection')(
            act)

    if gc.zero_init:
      init_w = hk.initializers.Constant(0.0)
    else:
      init_w = hk.initializers.VarianceScaling(scale=2., mode='fan_in')

    output_w = hk.get_parameter(
        'output_w',
        shape=(c.num_outer_channel, c.num_outer_channel,
               self.num_output_channel),
        dtype=act.dtype,
        init=init_w)
    output_b = hk.get_parameter(
        'output_b', shape=(self.num_output_channel,),
        dtype=act.dtype,
        init=hk.initializers.Constant(0.0))

    def compute_chunk(left_act):
      # This is equivalent to
      #
      # act = jnp.einsum('abc,ade->dceb', left_act, right_act)
      # act = jnp.einsum('dceb,cef->bdf', act, output_w) + output_b
      #
      # but faster.
      left_act = jnp.transpose(left_act, [0, 2, 1])
      act = jnp.einsum('acb,ade->dceb', left_act, right_act)
      act = jnp.einsum('dceb,cef->dbf', act, output_w) + output_b
      return jnp.transpose(act, [1, 0, 2])

    act = mapping.inference_subbatch(
        compute_chunk,
        c.chunk_size,
        batched_args=[left_act],
        nonbatched_args=[],
        low_memory=True,
        input_subbatch_dim=1,
        output_subbatch_dim=0)

    epsilon = 1e-3
    norm = jnp.einsum('abc,adc->bdc', mask, mask)
    act /= epsilon + norm

    return act

def dgram_from_positions(positions, num_bins, min_bin, max_bin):
  """Compute distogram from amino acid positions.
  Arguments:
    positions: [N_res, 3] Position coordinates.
    num_bins: The number of bins in the distogram.
    min_bin: The left edge of the first bin.
    max_bin: The left edge of the final bin. The final bin catches
        everything larger than `max_bin`.
  Returns:
    Distogram with the specified number of bins.
  """
  def squared_difference(x, y):
    return jnp.square(x - y)

  lower_breaks = jnp.linspace(min_bin, max_bin, num_bins)
  lower_breaks = jnp.square(lower_breaks)
  upper_breaks = jnp.concatenate([lower_breaks[1:],jnp.array([1e8], dtype=jnp.float32)], axis=-1)
  dist2 = jnp.sum(
      squared_difference(
          jnp.expand_dims(positions, axis=-2),
          jnp.expand_dims(positions, axis=-3)),
      axis=-1, keepdims=True)

  return ((dist2 > lower_breaks).astype(jnp.float32) * (dist2 < upper_breaks).astype(jnp.float32))

def dgram_from_positions_soft(positions, num_bins, min_bin, max_bin, temp=2.0):
  '''soft positions to dgram converter'''
  lower_breaks = jnp.append(-1e8,jnp.linspace(min_bin, max_bin, num_bins))
  upper_breaks = jnp.append(lower_breaks[1:],1e8)
  dist = jnp.sqrt(jnp.square(positions[...,:,None,:] - positions[...,None,:,:]).sum(-1,keepdims=True) + 1e-8)
  o = jax.nn.sigmoid((dist - lower_breaks)/temp) * jax.nn.sigmoid((upper_breaks - dist)/temp)
  o = o/(o.sum(-1,keepdims=True) + 1e-8)
  return o[...,1:]

def pseudo_beta_fn(aatype, all_atom_positions, all_atom_mask):
  """Create pseudo beta features."""
  
  ca_idx = residue_constants.atom_order['CA']
  cb_idx = residue_constants.atom_order['CB']

  is_gly = jnp.equal(aatype, residue_constants.restype_order['G'])
  is_gly_tile = jnp.tile(is_gly[..., None], [1] * len(is_gly.shape) + [3])
  pseudo_beta = jnp.where(is_gly_tile, all_atom_positions[..., ca_idx, :], all_atom_positions[..., cb_idx, :])

  if all_atom_mask is not None:
    pseudo_beta_mask = jnp.where(is_gly, all_atom_mask[..., ca_idx], all_atom_mask[..., cb_idx])
    pseudo_beta_mask = pseudo_beta_mask.astype(jnp.float32)
    return pseudo_beta, pseudo_beta_mask
  else:
    return pseudo_beta

class EvoformerIteration(hk.Module):
  """Single iteration (block) of Evoformer stack.
  Jumper et al. (2021) Suppl. Alg. 6 "EvoformerStack" lines 2-10
  """

  def __init__(self, config, global_config, is_extra_msa,
               name='evoformer_iteration'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config
    self.is_extra_msa = is_extra_msa

  def __call__(self, activations, masks, use_dropout, safe_key=None):
    """Builds EvoformerIteration module.

    Arguments:
      activations: Dictionary containing activations:
        * 'msa': MSA activations, shape [N_seq, N_res, c_m].
        * 'pair': pair activations, shape [N_res, N_res, c_z].
      masks: Dictionary of masks:
        * 'msa': MSA mask, shape [N_seq, N_res].
        * 'pair': pair mask, shape [N_res, N_res].
      safe_key: prng.SafeKey encapsulating rng key.

    Returns:
      Outputs, same shape/type as act.
    """
    c = self.config
    gc = self.global_config

    msa_act, pair_act = activations['msa'], activations['pair']

    if safe_key is None:
      safe_key = prng.SafeKey(hk.next_rng_key())

    msa_mask, pair_mask = masks['msa'], masks['pair']

    dropout_wrapper_fn = functools.partial(
        dropout_wrapper,
        global_config=gc,
        use_dropout=use_dropout)

    safe_key, *sub_keys = safe_key.split(10)
    sub_keys = iter(sub_keys)

    msa_act = dropout_wrapper_fn(
        MSARowAttentionWithPairBias(
            c.msa_row_attention_with_pair_bias, gc,
            name='msa_row_attention_with_pair_bias'),
        msa_act,
        msa_mask,
        safe_key=next(sub_keys),
        pair_act=pair_act)

    if not self.is_extra_msa:
      attn_mod = MSAColumnAttention(
          c.msa_column_attention, gc, name='msa_column_attention')
    else:
      attn_mod = MSAColumnGlobalAttention(
          c.msa_column_attention, gc, name='msa_column_global_attention')
    msa_act = dropout_wrapper_fn(
        attn_mod,
        msa_act,
        msa_mask,
        safe_key=next(sub_keys))

    msa_act = dropout_wrapper_fn(
        Transition(c.msa_transition, gc, name='msa_transition'),
        msa_act,
        msa_mask,
        safe_key=next(sub_keys))

    pair_act = dropout_wrapper_fn(
        OuterProductMean(
            config=c.outer_product_mean,
            global_config=self.global_config,
            num_output_channel=int(pair_act.shape[-1]),
            name='outer_product_mean'),
        msa_act,
        msa_mask,
        safe_key=next(sub_keys),
        output_act=pair_act)

    pair_act = dropout_wrapper_fn(
        TriangleMultiplication(c.triangle_multiplication_outgoing, gc,
                               name='triangle_multiplication_outgoing'),
        pair_act,
        pair_mask,
        safe_key=next(sub_keys))
    pair_act = dropout_wrapper_fn(
        TriangleMultiplication(c.triangle_multiplication_incoming, gc,
                               name='triangle_multiplication_incoming'),
        pair_act,
        pair_mask,
        safe_key=next(sub_keys))

    pair_act = dropout_wrapper_fn(
        TriangleAttention(c.triangle_attention_starting_node, gc,
                          name='triangle_attention_starting_node'),
        pair_act,
        pair_mask,
        safe_key=next(sub_keys))
    pair_act = dropout_wrapper_fn(
        TriangleAttention(c.triangle_attention_ending_node, gc,
                          name='triangle_attention_ending_node'),
        pair_act,
        pair_mask,
        safe_key=next(sub_keys))

    pair_act = dropout_wrapper_fn(
        Transition(c.pair_transition, gc, name='pair_transition'),
        pair_act,
        pair_mask,
        safe_key=next(sub_keys))

    return {'msa': msa_act, 'pair': pair_act}


class EmbeddingsAndEvoformer(hk.Module):
  """Embeds the input data and runs Evoformer.

  Produces the MSA, single and pair representations.
  Jumper et al. (2021) Suppl. Alg. 2 "Inference" line 5-18
  """

  def __init__(self, config, global_config, name='evoformer'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, batch, safe_key=None):

    c = self.config
    gc = self.global_config
    dtype = jnp.bfloat16 if gc.bfloat16 else jnp.float32

    if safe_key is None:
      safe_key = prng.SafeKey(hk.next_rng_key())
    
    with utils.bfloat16_context():

      # Embed clustered MSA.
      # Jumper et al. (2021) Suppl. Alg. 2 "Inference" line 5
      # Jumper et al. (2021) Suppl. Alg. 3 "InputEmbedder"

      msa_feat = batch['msa_feat'].astype(dtype)
      target_feat = jnp.pad(batch["target_feat"].astype(dtype),[[0,0],[1,1]])
      preprocess_1d = common_modules.Linear(c.msa_channel, name='preprocess_1d')(target_feat)
      preprocess_msa = common_modules.Linear(c.msa_channel, name='preprocess_msa')(msa_feat)
      msa_activations = preprocess_1d[None] + preprocess_msa

      left_single = common_modules.Linear(c.pair_channel, name='left_single')(target_feat)
      right_single = common_modules.Linear(c.pair_channel, name='right_single')(target_feat)
      pair_activations = left_single[:, None] + right_single[None]

      mask_2d = batch['seq_mask'][:, None] * batch['seq_mask'][None, :]
      mask_2d = mask_2d.astype(dtype)

      # Inject previous outputs for recycling.
      # Jumper et al. (2021) Suppl. Alg. 2 "Inference" line 6
      # Jumper et al. (2021) Suppl. Alg. 32 "RecyclingEmbedder"

      if gc.use_dgram:
        # use predicted distogram input (from Sergey)
        dgram = jax.nn.softmax(batch["prev_dgram"])
        dgram_map = jax.nn.one_hot(jnp.repeat(jnp.append(0,jnp.arange(15)),4),15).at[:,0].set(0)
        dgram = dgram @ dgram_map

      else:
        # use predicted position input
        prev_pseudo_beta = pseudo_beta_fn(batch['aatype'], batch['prev_pos'], None)
        if c.backprop_dgram:
          dgram = dgram_from_positions_soft(prev_pseudo_beta, temp=c.backprop_dgram_temp, **c.prev_pos)
        else:
          dgram = dgram_from_positions(prev_pseudo_beta, **c.prev_pos)    
      dgram = dgram.astype(dtype)    
      pair_activations += common_modules.Linear(c.pair_channel, name='prev_pos_linear')(dgram)

      if c.recycle_features:
        if 'prev_msa_first_row' in batch:
          prev_msa_first_row = common_modules.LayerNorm(
            axis=[-1], create_scale=True, create_offset=True,
            name='prev_msa_first_row_norm')(batch['prev_msa_first_row']).astype(dtype)
          msa_activations = msa_activations.at[0].add(prev_msa_first_row)

        if 'prev_pair' in batch:
          pair_activations += common_modules.LayerNorm(
            axis=[-1], create_scale=True, create_offset=True,
            name='prev_pair_norm')(batch['prev_pair']).astype(dtype)

      # Relative position encoding.
      # Jumper et al. (2021) Suppl. Alg. 4 "relpos"
      # Jumper et al. (2021) Suppl. Alg. 5 "one_hot"
      if c.max_relative_feature:
        # Add one-hot-encoded clipped residue distances to the pair activations.
        if "rel_pos" in batch:
          rel_pos = batch['rel_pos'].astype(dtype)
        else:
          if "offset" in batch:
            offset = batch['offset']
          else:
            pos = batch['residue_index']
            offset = pos[:, None] - pos[None, :]
          rel_pos = jax.nn.one_hot(
              jnp.clip(
                  offset + c.max_relative_feature,
                  a_min=0,
                  a_max=2 * c.max_relative_feature),
              2 * c.max_relative_feature + 1).astype(dtype)
        pair_activations += common_modules.Linear(c.pair_channel, name='pair_activiations')(rel_pos)

      # Embed templates into the pair activations.
      # Jumper et al. (2021) Suppl. Alg. 2 "Inference" lines 9-13

      if c.template.enabled:
        template_batch = {k: batch[k] for k in batch if k.startswith('template_')}

        multichain_mask = batch['asym_id'][:, None] == batch['asym_id'][None, :]
        multichain_mask = jnp.where(batch["mask_template_interchain"], multichain_mask, True)

        template_pair_representation = TemplateEmbedding(c.template, gc)(
            pair_activations,
            template_batch,
            mask_2d,
            multichain_mask,
            use_dropout=batch["use_dropout"])

        pair_activations += template_pair_representation

      # Embed extra MSA features.
      # Jumper et al. (2021) Suppl. Alg. 2 "Inference" lines 14-16
      if c.use_extra_msa:
        extra_msa_feat = create_extra_msa_feature(batch)
        extra_msa_activations = common_modules.Linear(c.extra_msa_channel,
                                                      name='extra_msa_activations')(extra_msa_feat).astype(dtype)
        # Extra MSA Stack.
        # Jumper et al. (2021) Suppl. Alg. 18 "ExtraMsaStack"
        extra_msa_stack_input = {'msa': extra_msa_activations,
                                 'pair': pair_activations}
        extra_msa_stack_iteration = EvoformerIteration(c.evoformer, gc,
                                                       is_extra_msa=True, name='extra_msa_stack')
        def extra_msa_stack_fn(x):
          act, safe_key = x
          safe_key, safe_subkey = safe_key.split()
          extra_evoformer_output = extra_msa_stack_iteration(
              activations=act,
              masks={'msa': batch['extra_msa_mask'].astype(dtype),
                     'pair': mask_2d},
              safe_key=safe_subkey,
              use_dropout=batch["use_dropout"])
          return (extra_evoformer_output, safe_key)
        if gc.use_remat: extra_msa_stack_fn = hk.remat(extra_msa_stack_fn)
        extra_msa_stack = layer_stack.layer_stack(c.extra_msa_stack_num_block)(extra_msa_stack_fn)
        extra_msa_output, safe_key = extra_msa_stack((extra_msa_stack_input, safe_key))
        pair_activations = extra_msa_output['pair']

      evoformer_input = {'msa': msa_activations,'pair': pair_activations}
      evoformer_masks = {'msa': batch['msa_mask'].astype(dtype),
                         'pair': mask_2d}
      ####################################################################

      # Append num_templ rows to msa_activations with template embeddings.
      # Jumper et al. (2021) Suppl. Alg. 2 "Inference" lines 7-8
      if c.template.enabled and c.template.embed_torsion_angles:
        num_templ, num_res = batch['template_aatype'].shape
        # Embed the templates aatypes.
        aatype = batch['template_aatype']
        aatype_one_hot = jax.nn.one_hot(batch['template_aatype'], 22, axis=-1)

        # Embed the templates aatype, torsion angles and masks.
        # Shape (templates, residues, msa_channels)
        ret = all_atom.atom37_to_torsion_angles(
            aatype=aatype,
            all_atom_pos=batch['template_all_atom_positions'],
            all_atom_mask=batch['template_all_atom_mask'],
            # Ensure consistent behaviour during testing:
            placeholder_for_undefined=not gc.zero_init)

        template_features = jnp.concatenate([
            aatype_one_hot,
            jnp.reshape(ret['torsion_angles_sin_cos'], [num_templ, num_res, 14]),
            jnp.reshape(ret['alt_torsion_angles_sin_cos'], [num_templ, num_res, 14]),
            ret['torsion_angles_mask']], axis=-1).astype(dtype)

        template_activations = common_modules.Linear(
            c.msa_channel,
            initializer='relu',
            name='template_single_embedding')(template_features)
        template_activations = jax.nn.relu(template_activations)
        template_activations = common_modules.Linear(
            c.msa_channel,
            initializer='relu',
            name='template_projection')(template_activations)

        # Concatenate the templates to the msa.
        evoformer_input['msa'] = jnp.concatenate([evoformer_input['msa'], template_activations], axis=0)

        # Concatenate templates masks to the msa masks.
        # Use mask from the psi angle, as it only depends on the backbone atoms
        # from a single residue.
        torsion_angle_mask = ret['torsion_angles_mask'][:, :, 2]
        torsion_angle_mask = torsion_angle_mask.astype(evoformer_masks['msa'].dtype)
        evoformer_masks['msa'] = jnp.concatenate([evoformer_masks['msa'], torsion_angle_mask], axis=0)    
      ####################################################################
      if c.use_msa:
        # Main trunk of the network
        # Jumper et al. (2021) Suppl. Alg. 2 "Inference" lines 17-18
        evoformer_iteration = EvoformerIteration(c.evoformer, gc, is_extra_msa=False, name='evoformer_iteration')
        def evoformer_fn(x):
          act, safe_key = x
          safe_key, safe_subkey = safe_key.split()
          evoformer_output = evoformer_iteration(
              activations=act,
              masks=evoformer_masks,
              safe_key=safe_subkey,
              use_dropout=batch["use_dropout"])
          return (evoformer_output, safe_key)
        if gc.use_remat: evoformer_fn = hk.remat(evoformer_fn)
        evoformer_stack = layer_stack.layer_stack(c.evoformer_num_block)(evoformer_fn)
        evoformer_output, safe_key = evoformer_stack((evoformer_input, safe_key))
        msa_activations = evoformer_output['msa']
        pair_activations = evoformer_output['pair']

      single_activations = common_modules.Linear(c.seq_channel, name='single_activations')(msa_activations[0])
      num_sequences = batch['msa_feat'].shape[0]
      
    output = {
        'single': single_activations,
        'pair': pair_activations,
        # Crop away template rows such that they are not used in MaskedMsaHead.
        'msa': msa_activations[:num_sequences, :, :],
        'msa_first_row': msa_activations[0],
    }
    
    # Convert back to float32 if we're not saving memory.
    if not gc.bfloat16_output:
      for k, v in output.items():
        if v.dtype == jnp.bfloat16:
          output[k] = v.astype(jnp.float32)
    return output
  
####################################################################
####################################################################
class SingleTemplateEmbedding(hk.Module):
  """Embeds a single template.
  Jumper et al. (2021) Suppl. Alg. 2 "Inference" lines 9+11
  """

  def __init__(self, config, global_config, name='single_template_embedding'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, query_embedding, batch, mask_2d, multichain_mask_2d, use_dropout):
    """Build the single template embedding.
    Arguments:
      query_embedding: Query pair representation, shape [N_res, N_res, c_z].
      batch: A batch of template features (note the template dimension has been
        stripped out as this module only runs over a single template).
      mask_2d: Padding mask (Note: this doesn't care if a template exists,
        unlike the template_pseudo_beta_mask).
    Returns:
      A template embedding [N_res, N_res, c_z].
    """
    assert mask_2d.dtype == query_embedding.dtype
    dtype = query_embedding.dtype
    num_res = batch['template_aatype'].shape[0]
    num_channels = (self.config.template_pair_stack
                    .triangle_attention_ending_node.value_dim)

    if "template_dgram" in batch:
      template_dgram = batch["template_dgram"]
      template_mask_2d = batch["template_dgram"].sum(-1)

    else:
      template_mask = batch['template_pseudo_beta_mask']
      template_mask_2d = template_mask[:, None] * template_mask[None, :]
      if self.config.backprop_dgram:
        template_dgram = dgram_from_positions_soft(batch['template_pseudo_beta'],
                                                   temp=self.config.backprop_dgram_temp,
                                                   **self.config.dgram_features)
      else:
        template_dgram = dgram_from_positions(batch['template_pseudo_beta'],
                                              **self.config.dgram_features)

    template_mask_2d = template_mask_2d * multichain_mask_2d
    template_mask_2d = template_mask_2d.astype(dtype)

    template_dgram *= template_mask_2d[..., None]
    template_dgram = template_dgram.astype(dtype)
    to_concat = [template_dgram, template_mask_2d[:, :, None]]

    aatype = jax.nn.one_hot(batch['template_aatype'], 22, axis=-1, dtype=dtype)

    to_concat.append(jnp.tile(aatype[None, :, :], [num_res, 1, 1]))
    to_concat.append(jnp.tile(aatype[:, None, :], [1, num_res, 1]))

    if "template_dgram" in batch:
      unit_vector = [jnp.zeros((num_res,num_res,1))] * 3

    else:
      # Backbone affine mask: whether the residue has C, CA, N
      # (the template mask defined above only considers pseudo CB).
      n, ca, c = [residue_constants.atom_order[a] for a in ('N', 'CA', 'C')]
      template_mask = (
          batch['template_all_atom_mask'][..., n] *
          batch['template_all_atom_mask'][..., ca] *
          batch['template_all_atom_mask'][..., c])
      template_mask_2d = template_mask[:, None] * template_mask[None, :]
      template_mask_2d = template_mask_2d * multichain_mask_2d

      # compute unit_vector (not used by default)
      if self.config.use_template_unit_vector:
        raw_atom_pos = template_batch["template_all_atom_positions"]
        if gc.bfloat16:
          raw_atom_pos = raw_atom_pos.astype(jnp.float32)

        rot, trans = quat_affine.make_transform_from_reference(
            n_xyz=raw_atom_pos[:, n],
            ca_xyz=raw_atom_pos[:, ca],
            c_xyz=raw_atom_pos[:, c])
        affines = quat_affine.QuatAffine(
            quaternion=quat_affine.rot_to_quat(rot, unstack_inputs=True),
            translation=trans,
            rotation=rot,
            unstack_inputs=True)
        points = [jnp.expand_dims(x, axis=-2) for x in affines.translation]
        affine_vec = affines.invert_point(points, extra_dims=1)
        inv_distance_scalar = jax.lax.rsqrt(1e-6 + sum([jnp.square(x) for x in affine_vec]))
        inv_distance_scalar *= template_mask_2d.astype(inv_distance_scalar.dtype)
        unit_vector = [(x * inv_distance_scalar)[..., None] for x in affine_vec]
      else:      
        unit_vector = [jnp.zeros((num_res,num_res,1))] * 3
      
    unit_vector = [x.astype(dtype) for x in unit_vector]
    to_concat.extend(unit_vector)

    template_mask_2d = template_mask_2d.astype(dtype)
    to_concat.append(template_mask_2d[..., None])

    act = jnp.concatenate(to_concat, axis=-1)

    # Mask out non-template regions so we don't get arbitrary values in the
    # distogram for these regions.
    act *= template_mask_2d[..., None]

    # Jumper et al. (2021) Suppl. Alg. 2 "Inference" line 9
    act = common_modules.Linear(
        num_channels,
        initializer='relu',
        name='embedding2d')(act)

    # Jumper et al. (2021) Suppl. Alg. 2 "Inference" line 11
    act = TemplatePairStack(
        self.config.template_pair_stack, self.global_config)(act, mask_2d, use_dropout=use_dropout)

    act = common_modules.LayerNorm([-1], True, True, name='output_layer_norm')(act)
    return act


class TemplateEmbedding(hk.Module):
  """Embeds a set of templates.
  Jumper et al. (2021) Suppl. Alg. 2 "Inference" lines 9-12
  Jumper et al. (2021) Suppl. Alg. 17 "TemplatePointwiseAttention"
  """

  def __init__(self, config, global_config, name='template_embedding'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, query_embedding, template_batch, mask_2d, multichain_mask_2d, use_dropout):
    """Build TemplateEmbedding module.
    Arguments:
      query_embedding: Query pair representation, shape [N_res, N_res, c_z].
      template_batch: A batch of template features.
      mask_2d: Padding mask (Note: this doesn't care if a template exists,
        unlike the template_pseudo_beta_mask).
    Returns:
      A template embedding [N_res, N_res, c_z].
    """

    num_templates = template_batch['template_mask'].shape[0]
    num_channels = (self.config.template_pair_stack
                    .triangle_attention_ending_node.value_dim)
    num_res = query_embedding.shape[0]

    dtype = query_embedding.dtype
    template_mask = template_batch['template_mask']
    template_mask = template_mask.astype(dtype)

    query_num_channels = query_embedding.shape[-1]

    # Make sure the weights are shared across templates by constructing the
    # embedder here.
    # Jumper et al. (2021) Suppl. Alg. 2 "Inference" lines 9-12
    template_embedder = SingleTemplateEmbedding(self.config, self.global_config)

    def map_fn(batch):
      return template_embedder(query_embedding, batch, mask_2d, multichain_mask_2d,
        use_dropout=use_dropout)

    template_pair_representation = mapping.sharded_map(map_fn, in_axes=0)(template_batch)

    # Cross attend from the query to the templates along the residue
    # dimension by flattening everything else into the batch dimension.
    # Jumper et al. (2021) Suppl. Alg. 17 "TemplatePointwiseAttention"
    flat_query = jnp.reshape(query_embedding,[num_res * num_res, 1, query_num_channels])

    flat_templates = jnp.reshape(
        jnp.transpose(template_pair_representation, [1, 2, 0, 3]),
        [num_res * num_res, num_templates, num_channels])

    bias = (1e9 * (template_mask[None, None, None, :] - 1.))

    template_pointwise_attention_module = Attention(
        self.config.attention, self.global_config, query_num_channels)
    nonbatched_args = [bias]
    batched_args = [flat_query, flat_templates]

    embedding = mapping.inference_subbatch(
        template_pointwise_attention_module,
        self.config.subbatch_size,
        batched_args=batched_args,
        nonbatched_args=nonbatched_args,
        low_memory=self.config.subbatch_size is not None)
    embedding = jnp.reshape(embedding,[num_res, num_res, query_num_channels])

    # No gradients if no templates.
    embedding *= (jnp.sum(template_mask) > 0.).astype(embedding.dtype)

    return embedding
####################################################################