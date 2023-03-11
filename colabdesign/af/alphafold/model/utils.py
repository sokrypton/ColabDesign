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

"""A collection of JAX utility functions for use in protein folding."""

import collections
import contextlib
import functools
import numbers
from typing import Mapping

import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np
import io

def bfloat16_creator(next_creator, shape, dtype, init, context):
  """Creates float32 variables when bfloat16 is requested."""
  if context.original_dtype == jnp.bfloat16:
    dtype = jnp.float32
  return next_creator(shape, dtype, init)

def bfloat16_getter(next_getter, value, context):
  """Casts float32 to bfloat16 when bfloat16 was originally requested."""
  if context.original_dtype == jnp.bfloat16:
    assert value.dtype == jnp.float32
    value = value.astype(jnp.bfloat16)
  return next_getter(value)

@contextlib.contextmanager
def bfloat16_context():
  with hk.custom_creator(bfloat16_creator), hk.custom_getter(bfloat16_getter):
    yield

def final_init(config):
  if config.zero_init:
    return 'zeros'
  else:
    return 'linear'

def batched_gather(params, indices, axis=0, batch_dims=0):
  """Implements a JAX equivalent of `tf.gather` with `axis` and `batch_dims`."""
  take_fn = lambda p, i: jnp.take(p, i, axis=axis, mode="clip")
  for _ in range(batch_dims):
    take_fn = jax.vmap(take_fn)
  return take_fn(params, indices)


def mask_mean(mask, value, axis=None, drop_mask_channel=False, eps=1e-10):
  """Masked mean."""
  if drop_mask_channel:
    mask = mask[..., 0]

  mask_shape = mask.shape
  value_shape = value.shape

  assert len(mask_shape) == len(value_shape)

  if isinstance(axis, numbers.Integral):
    axis = [axis]
  elif axis is None:
    axis = list(range(len(mask_shape)))

  broadcast_factor = 1.
  for axis_ in axis:
    value_size = value_shape[axis_]
    mask_size = mask_shape[axis_]
    if mask_size == 1:
      broadcast_factor *= value_size
    else:
      assert mask_size == value_size

  return (jnp.sum(mask * value, axis=axis) /
          (jnp.sum(mask, axis=axis) * broadcast_factor + eps))

def flat_params_to_haiku(params, fuse=None):
  """Convert a dictionary of NumPy arrays to Haiku parameters."""
  P = {}
  for path, array in params.items():
    scope, name = path.split('//')
    if scope not in P:
      P[scope] = {}
    P[scope][name] = array
  if fuse is not None:
    for a in ["evoformer_iteration",
              "extra_msa_stack",
              "template_embedding/single_template_embedding/template_embedding_iteration",
              "template_embedding/single_template_embedding/template_pair_stack/__layer_stack_no_state"]:
      for b in ["triangle_multiplication_incoming","triangle_multiplication_outgoing"]:
        k = f"alphafold/alphafold_iteration/evoformer/{a}/{b}"

        if fuse and f"{k}/center_layer_norm" in P:
          for c in ["gate","projection"]:
            L = P.pop(f"{k}/left_{c}")
            R = P.pop(f"{k}/right_{c}")
            P[f"{k}/{c}"] = {}
            for d in ["bias","weights"]:
              P[f"{k}/{c}"][d] = np.concatenate([L[d],R[d]],-1)
          P[f"{k}/center_norm"] = P.pop(f"{k}/center_layer_norm")
          P[f"{k}/left_norm_input"] = P.pop(f"{k}/layer_norm_input")

        if not fuse and f"{k}/center_norm" in P:
          for c in ["gate","projection"]:
            LR = P.pop(f"{k}/{c}")
            P[f"{k}/left_{c}"] = {}
            P[f"{k}/right_{c}"] = {}
            for d in ["bias","weights"]:
              half = LR[d].shape[-1] // 2
              P[f"{k}/left_{c}"][d] = LR[d][...,:half]
              P[f"{k}/right_{c}"][d] = LR[d][...,half:]
          P[f"{k}/center_layer_norm"] = P.pop(f"{k}/center_norm")
          P[f"{k}/layer_norm_input"] = P.pop(f"{k}/left_norm_input")
  return P