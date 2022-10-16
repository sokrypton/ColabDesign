# Copyright (c) Levinthal, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import jax
import math
import jax.numpy as jnp
import haiku as hk
from jax.nn import softmax

from colabdesign.shared.prng import SafeKey

class RowSelfAttention(hk.Module):
  """Compute self-attention over rows of a 2D input."""

  def __init__(
    self,
    config,
  ):
    super().__init__()
    self.head_num = config.RowAtt.head_num
    self.embed_dim = config.RowAtt.embed_dim
    self.dropout = config.dropout
    self.max_tokens_per_msa = config.max_tokens_per_msa

    self.head_dim = self.embed_dim // self.head_num
    self.scaling = self.head_dim ** -0.5
    self.attn_shape = "hij"

    self.k_proj = hk.Linear(self.embed_dim, name="k_proj")
    self.v_proj = hk.Linear(self.embed_dim, name="v_proj")
    self.q_proj = hk.Linear(self.embed_dim, name="q_proj")
    self.out_proj = hk.Linear(self.embed_dim, name="out_proj")
    self.safe_key = SafeKey(hk.next_rng_key())

  def align_scaling(self, q):
    num_rows = q.shape[0]
    return self.scaling / math.sqrt(num_rows)

  def _batched_forward(
    self,
    x,
    self_attn_padding_mask,
  ):
    num_rows, num_cols, embed_dim = x.shape
    max_rows = max(1, self.max_tokens_per_msa // num_cols)
    attns = 0
    scaling = self.align_scaling(x)
    for start in range(0, num_rows, max_rows):
      attn_weights = self.compute_attention_weights(
        x[start: start + max_rows],
        scaling,
        self_attn_padding_mask=self_attn_padding_mask[:, start: start + max_rows]
      )
      attns += attn_weights
    attn_probs = softmax(attns, -1)
    self.safe_key, use_key = self.safe_key.split()
    attn_probs = hk.dropout(use_key.get(), self.dropout, attn_probs)

    outputs = []
    for start in range(0, num_rows, max_rows):
      output = self.compute_attention_update(x[start: start + max_rows], attn_probs)
      outputs.append(output)

    output = jnp.concatenate(outputs, 0)
    return output, attn_probs

  def compute_attention_weights(
    self,
    x,
    scaling: float,
    self_attn_padding_mask,
  ):
    num_rows, num_cols, embed_dim = x.shape
    q = self.q_proj(x).reshape([num_rows, num_cols, self.head_num, self.head_dim])
    k = self.k_proj(x).reshape([num_rows, num_cols, self.head_num, self.head_dim])
    q *= scaling
    # Zero out any padded aligned positions - this is important since
    # we take a sum across the alignment axis.
    q *= 1 - jnp.expand_dims(jnp.expand_dims(self_attn_padding_mask, 2), 3)

    attn_weights = jnp.einsum(f"rihd,rjhd->{self.attn_shape}", q, k)
    attn_weights *= 1 - jnp.expand_dims(jnp.expand_dims(self_attn_padding_mask[0], 0), 2)
    attn_weights += jnp.expand_dims(jnp.expand_dims(self_attn_padding_mask[0], 0), 2) * -10000

    return attn_weights

  def compute_attention_update(
    self,
    x,
    attn_probs,
  ):
    num_rows, num_cols, embed_dim = x.shape
    v = self.v_proj(x).reshape([num_rows, num_cols, self.head_num, self.head_dim])
    context = jnp.einsum(f"{self.attn_shape},rjhd->rihd", attn_probs, v)
    context = context.reshape([num_rows, num_cols, embed_dim])
    output = self.out_proj(context)
    return output

  def __call__(self, x,
         self_attn_padding_mask,):
    
    num_rows, num_cols, embed_dim = x.shape
    if num_rows * num_cols > self.max_tokens_per_msa:
      return self._batched_forward(x, self_attn_padding_mask)
    else:
      scaling = self.align_scaling(x)
      attn_weights = self.compute_attention_weights(
        x, scaling, self_attn_padding_mask
      )
      attn_probs = softmax(attn_weights, -1)
      self.safe_key, use_key = self.safe_key.split()
      attn_probs = hk.dropout(use_key.get(), self.dropout, attn_probs)
      output = self.compute_attention_update(x, attn_probs)
      return output, attn_probs


class ColumnSelfAttention(hk.Module):
  """Compute self-attention over columns of a 2D input."""

  def __init__(self, config):
    super().__init__()
    self.head_num = config.ColAtt.head_num
    self.embed_dim = config.RowAtt.embed_dim
    self.dropout = config.dropout
    self.max_tokens_per_msa = config.max_tokens_per_msa

    self.head_dim = self.embed_dim // self.head_num
    self.scaling = self.head_dim ** -0.5
    self.safe_key = SafeKey(hk.next_rng_key())

    self.k_proj = hk.Linear(self.embed_dim, name='k_proj')
    self.v_proj = hk.Linear(self.embed_dim, name='v_proj')
    self.q_proj = hk.Linear(self.embed_dim, name='q_proj')
    self.out_proj = hk.Linear(self.embed_dim, name='out_proj')

  def _batched_forward(
    self,
    x,
    self_attn_padding_mask,
  ):
    num_rows, num_cols, embed_dim = x.shape
    max_cols = max(1, self.max_tokens_per_msa // num_rows)
    outputs = []
    attns = []
    for start in range(0, num_cols, max_cols):
      output, attn = self.compute_attention_update(
        x[:, start: start + max_cols],
        self_attn_padding_mask=self_attn_padding_mask[:, :, start: start + max_cols]
      )
      outputs.append(output)
      attns.append(attn)
    output = jnp.concatenate(outputs, 1)
    attns = jnp.concatenate(attns, 1)
    return output, attns

  def compute_attention_update(
    self,
    x,
    self_attn_padding_mask,
  ):
    num_rows, num_cols, embed_dim = x.shape
    if num_rows == 1:
      attn_probs = jnp.ones(
        [self.head_num, num_cols, num_rows, num_rows],
        dtype=x.dtype,
      )
      output = self.out_proj(self.v_proj(x))
      return output, attn_probs
    else:
      q = self.q_proj(x).reshape([num_rows, num_cols, self.head_num, self.head_dim])
      k = self.k_proj(x).reshape([num_rows, num_cols, self.head_num, self.head_dim])
      v = self.v_proj(x).reshape([num_rows, num_cols, self.head_num, self.head_dim])
      q *= self.scaling

      attn_weights = jnp.einsum("ichd,jchd->hcij", q, k)
      attn_weights *= 1 - jnp.expand_dims(jnp.expand_dims(self_attn_padding_mask.transpose(), 0), 3)
      attn_weights += jnp.expand_dims(jnp.expand_dims(self_attn_padding_mask.transpose(), 0), 3) * -10000
      attn_probs = softmax(attn_weights, -1)

      self.safe_key, use_key = self.safe_key.split()
      attn_probs = hk.dropout(use_key.get(), self.dropout, attn_probs)
      context = jnp.einsum("hcij,jchd->ichd", attn_probs, v)
      context = context.reshape([num_rows, num_cols, embed_dim])
      output = self.out_proj(context)
      return output, attn_probs

  def __call__(
    self,
    x,
    self_attn_padding_mask,
  ):
    # if False and num_rows * num_cols > 2 ** 14 and not torch.is_grad_enabled():
    num_rows, num_cols, embed_dim = x.shape
    if num_rows * num_cols > self.max_tokens_per_msa:
      return self._batched_forward(x, self_attn_padding_mask)
    else:
      return self.compute_attention_update(x, self_attn_padding_mask)

