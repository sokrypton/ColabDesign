# Copyright (c) Facebook, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import joblib
import jax.numpy as jnp
import numpy as np
import haiku as hk
import jax

from .modules import (
  AxialTransformerLayer,
  EmbedPosition,
  MSAPositionEmbedding,
  ContactPredictionHead,
  LmHead,
)

from colabdesign.shared.prng import SafeKey

class MSATransformer(hk.Module):
  def __init__(self, alphabet, config):
    super().__init__()
    self.alphabet_size = len(alphabet)
    self.padding_idx = alphabet.padding_idx
    self.mask_idx = alphabet.mask_idx
    self.cls_idx = alphabet.cls_idx
    self.eos_idx = alphabet.eos_idx
    self.prepend_bos = alphabet.prepend_bos
    self.append_eos = alphabet.append_eos
    self.config = config
    self.dropout = config.dropout

    self.embed_tokens = hk.Embed(
      vocab_size=self.alphabet_size,
      embed_dim=self.config.embed_dim, 
    )

    self.msa_position_embedding = MSAPositionEmbedding(self.config.embed_dim)
    self.safe_key = SafeKey(hk.next_rng_key())

    self.layers = [
      AxialTransformerLayer(self.config)
      for _ in range(self.config.layer_num)
    ]

    self.contact_head = ContactPredictionHead(
      self.config.layer_num * self.config.RowAtt.head_num,
      self.prepend_bos,
      self.append_eos,
      eos_idx=self.eos_idx,
    )
    self.embed_positions = EmbedPosition(
      self.config,
      self.padding_idx,
    )

    self.emb_layer_norm_before = hk.LayerNorm(-1, create_scale=True, create_offset=True)
    self.emb_layer_norm_after = hk.LayerNorm(-1, create_scale=True, create_offset=True)

    self.lm_head = LmHead(
      config=self.config,
      output_dim=self.alphabet_size,
      weight=self.embed_tokens.embeddings.transpose(),
    )

  def __call__(self, tokens):
    num_alignments, seqlen = tokens.shape
    padding_mask = jnp.equal(tokens, self.padding_idx)  # R, C
    x = self.embed_tokens(tokens)
    x += self.embed_positions(tokens)
    x += self.msa_position_embedding(tokens)
    x = self.emb_layer_norm_before(x)

    self.safe_key, use_key = self.safe_key.split()
    x = hk.dropout(use_key.get(), self.dropout, x)
    x = x * (1 - jnp.expand_dims(padding_mask, axis=-1))

    row_attn_weights = []
    col_attn_weights = []

    for layer_idx, layer in enumerate(self.layers):
      x = layer(
        x,
        self_attn_padding_mask=padding_mask,
      )
      x, col_attn, row_attn = x
      col_attn_weights.append(col_attn)
      row_attn_weights.append(row_attn)

    x = self.emb_layer_norm_after(x)
    x = self.lm_head(x)

    result = {"logits": x}
    # col_attentions: L x H x C x R x R
    col_attentions = jnp.stack(col_attn_weights, 0)
    # row_attentions: L x H x C x C
    row_attentions = jnp.stack(row_attn_weights, 0)
    result["col_attentions"] = col_attentions
    result["row_attentions"] = row_attentions
    contacts = self.contact_head(tokens, row_attentions)
    result["contacts"] = contacts

    return result


class RunModel:
  '''container for msa transformer'''

  def __init__(self, alphabet, config):
    self.padding_idx = alphabet.padding_idx

    def _forward(tokens):
      model = MSATransformer(alphabet, config)
      return model(tokens)

    _forward_t = hk.transform(_forward)
    self.init = jax.jit(_forward_t.init)
    self.apply = jax.jit(_forward_t.apply)
    self.key = jax.random.PRNGKey(42)

  def load_params(self, path):
    self.params = joblib.load(path)

  def __call__(self, tokens):
    assert tokens.ndim == 2
    num_alignments, seqlen = tokens.shape

    if num_alignments > 1024:
      raise RuntimeError(
        "Using model with MSA position embedding trained on maximum MSA "
        f"depth of 1024, but received {num_alignments} alignments."
      )

    self.key, use_key = jax.random.split(self.key)
    result = self.apply(self.params, use_key, tokens)
    result_new = {}
    for ikey in result.keys():
      result_new[ikey] = np.array(result[ikey])
    return result_new
