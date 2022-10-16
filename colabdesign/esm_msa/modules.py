# Copyright (c) Facebook, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from typing import Optional

import haiku as hk
import jax
import jax.numpy as jnp

from .axial_attention import ColumnSelfAttention, RowSelfAttention
from colabdesign.shared.prng import SafeKey



def symmetrize(x):
  "Make layer symmetric in final two dimensions, used for contact prediction."
  return x + x.transpose([0, 2, 1])


def apc(x):
  "Perform average product correct, used for contact prediction."
  a1 = x.sum(-1, keepdims=True)
  a2 = x.sum(-2, keepdims=True)
  a12 = x.sum((-1, -2), keepdims=True)

  avg = a1 * a2
  avg = avg / a12
  normalized = x - avg
  return normalized


class AxialTransformerLayer(hk.Module):
  """Implements an Axial MSA Transformer block."""

  def __init__(
    self,
    config,
  ) -> None:
    super().__init__()
    self.config = config

    row_self_attention = RowSelfAttention(config)
    column_self_attention = ColumnSelfAttention(config)
    feed_forward_layer = FeedForwardNetwork(config)

    self.row_self_attention = self.build_residual(row_self_attention, name='row_self_attention')
    self.column_self_attention = self.build_residual(column_self_attention, name='column_self_attention')
    self.feed_forward_layer = self.build_residual(feed_forward_layer, name='feed_forward_layer')

  def build_residual(self, layer: hk.Module, name=None):
    return NormalizedResidualBlock(
      layer,
      self.config,
      name=name,
    )

  def __call__(
    self,
    x,
    self_attn_padding_mask,
  ):
    """
    LayerNorm is applied either before or after the self-attention/ffn
    modules similar to the original Transformer implementation.
    """
    x, row_attn = self.row_self_attention(
      x,
      self_attn_padding_mask=self_attn_padding_mask,
    )
    x, column_attn = self.column_self_attention(
      x,
      self_attn_padding_mask=self_attn_padding_mask,
    )
    x = self.feed_forward_layer(x)
    return x, column_attn, row_attn


class LmHead(hk.Module):
  def __init__(self, config, output_dim, weight):
    super().__init__()
    self.layer_norm = hk.LayerNorm(-1, create_scale=True, create_offset=True)
    self.dense = hk.Linear(config.embed_dim, name='dense')
    self.weight = weight
    self.bias = hk.get_parameter(name='bias', shape=[output_dim], init=jnp.zeros)

  def __call__(self, input):
    x = self.dense(input)
    x = jax.nn.gelu(x)
    x = self.layer_norm(x)
    x = jnp.dot(x, self.weight) + self.bias
    return x


class ContactPredictionHead(hk.Module):
  """Performs symmetrization, apc, and computes a logistic regression on the output features"""

  def __init__(
    self,
    in_features: int,
    prepend_bos: bool,
    append_eos: bool,
    bias=True,
    eos_idx: Optional[int] = None,
  ):
    super().__init__()
    self.in_features = in_features
    self.prepend_bos = prepend_bos
    self.append_eos = append_eos
    self.eos_idx = eos_idx
    self.regression = hk.Linear(1, with_bias=bias)
    self.activation = jax.nn.sigmoid

  def __call__(self, tokens, attentions):
    # remove eos token attentions
    if self.append_eos:
      eos_mask = jnp.not_equal(tokens, self.eos_idx)
      eos_mask = jnp.expand_dims(eos_mask, axis=0) * jnp.expand_dims(eos_mask, axis=1)
      attentions = attentions * eos_mask[None, None, :, :]
      attentions = attentions[..., :-1, :-1]

    # remove cls token attentions
    if self.prepend_bos:
      attentions = attentions[..., 1:, 1:]

    layers, heads, seqlen, _ = attentions.shape
    attentions = attentions.reshape([layers * heads, seqlen, seqlen])

    # features: C x T x T
    attentions = apc(symmetrize(attentions))
    attentions = attentions.transpose([1, 2, 0])
    return self.activation(self.regression(attentions).squeeze(2))


class NormalizedResidualBlock(hk.Module):
  def __init__(
    self,
    layer: hk.Module,
    config,
    name=None,
  ):
    super().__init__(name=name)
    self.embed_dim = config.embed_dim
    self.dropout = config.dropout
    self.safe_key = SafeKey(hk.next_rng_key())

    self.layer = layer
    self.layer_norm = hk.LayerNorm(-1, create_scale=True, create_offset=True)

  def __call__(self, x, *args, **kwargs):
    residual = x
    x = self.layer_norm(x)
    outputs = self.layer(x, *args, **kwargs)
    if isinstance(outputs, tuple):
      x, *out = outputs
    else:
      x = outputs
      out = None

    self.safe_key, use_key = self.safe_key.split()
    x = hk.dropout(use_key.get(), self.dropout, x)
    x = residual + x

    if out is not None:
      return (x,) + tuple(out)
    else:
      return x


class FeedForwardNetwork(hk.Module):
  def __init__(
    self,
    config,
  ):
    super().__init__()
    self.embed_dim = config.embed_dim
    self.ffn_embed_dim = config.Ffn.embed_dim
    self.max_tokens_per_msa = config.max_tokens_per_msa
    self.dropout = config.dropout

    self.safe_key = SafeKey(hk.next_rng_key())
    self.activation_fn = jax.nn.gelu

    self.fc1 = hk.Linear(self.ffn_embed_dim, name='fc1')
    self.fc2 = hk.Linear(self.embed_dim, name='fc2')

  def __call__(self, x):
    x = self.activation_fn(self.fc1(x))
    self.safe_key, use_key = self.safe_key.split()
    x = hk.dropout(use_key.get(), self.dropout, x)
    x = self.fc2(x)
    return x


class MSAPositionEmbedding(hk.Module):
  def __init__(self, embed_dim):
    super().__init__()
    self.embed_dim = embed_dim
    self.weight = hk.get_parameter(name='data',
                     shape=[1024, 1, embed_dim],
                     init=jnp.zeros)

  def __call__(self, x):
    # num_alignments, seq_len = x.shape
    num_rows, num_cols = x.shape
    return self.weight[:num_rows]


class EmbedPosition(hk.Module):
  def __init__(self, config, padding_idx):
    super().__init__()
    self.max_position = config.max_position
    self.embed_dim = config.embed_dim
    self.padding_idx = padding_idx
    self.max_position_ = self.max_position + self.padding_idx + 1
    self.embed = hk.Embed(vocab_size=self.max_position_,
                embed_dim=self.embed_dim)

  def __call__(self, tokens):
    mask = jnp.not_equal(tokens, self.padding_idx)
    # tokens always begin with <cls>, do not consider. <cls> is before <pad> in alphabet.
    positions = jnp.cumsum(mask, axis=-1, dtype='int32') * mask + self.padding_idx

    # position_oh = jax.nn.one_hot(positions, self.max_position_)
    # weight = hk.get_parameter('weight', shape=[self.max_position_, self.embed_dim], init=jnp.zeros)
    # x = jnp.dot(position_one_hot, weight)
    # return x
    return self.embed(positions)
