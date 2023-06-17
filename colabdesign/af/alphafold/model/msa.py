###################
# TODO: msa.py
###################
#
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

"""
MSA sampling pipeline is moved inside the JAX model 
for easier implementation of recycling and ensembling.
"""

import functools
from typing import Sequence

import jax
import jax.numpy as jnp
from colabdesign.af.alphafold.model import utils

def gumbel_noise(key: jnp.ndarray, shape: Sequence[int]) -> jnp.ndarray:
  """Generate Gumbel Noise of given Shape.

  This generates samples from Gumbel(0, 1).

  Args:
    key: Jax random number key.
    shape: Shape of noise to return.

  Returns:
    Gumbel noise of given shape.
  """
  epsilon = 1e-6
  uniform = utils.padding_consistent_rng(jax.random.uniform)
  uniform_noise = uniform(
      key, shape=shape, dtype=jnp.float32, minval=0., maxval=1.)
  gumbel = -jnp.log(-jnp.log(uniform_noise + epsilon) + epsilon)
  return gumbel


def gumbel_max_sample(key: jnp.ndarray, logits: jnp.ndarray) -> jnp.ndarray:
  """Samples from a probability distribution given by 'logits'.

  This uses Gumbel-max trick to implement the sampling in an efficient manner.

  Args:
    key: prng key.
    logits: Logarithm of probabilities to sample from, probabilities can be
      unnormalized.

  Returns:
    Sample from logprobs in one-hot form.
  """
  z = gumbel_noise(key, logits.shape)
  return jax.nn.one_hot(
      jnp.argmax(logits + z, axis=-1),
      logits.shape[-1],
      dtype=logits.dtype)

def gumbel_argsort_sample_idx(key: jnp.ndarray,
                              logits: jnp.ndarray) -> jnp.ndarray:
  """Samples with replacement from a distribution given by 'logits'.

  This uses Gumbel trick to implement the sampling an efficient manner. For a
  distribution over k items this samples k times without replacement, so this
  is effectively sampling a random permutation with probabilities over the
  permutations derived from the logprobs.

  Args:
    key: prng key.
    logits: Logarithm of probabilities to sample from, probabilities can be
      unnormalized.

  Returns:
    Sample from logprobs in one-hot form.
  """
  z = gumbel_noise(key, logits.shape)
  # This construction is equivalent to jnp.argsort, but using a non stable sort,
  # since stable sort's aren't supported by jax2tf.
  axis = len(logits.shape) - 1
  iota = jax.lax.broadcasted_iota(jnp.int64, logits.shape, axis)
  _, perm = jax.lax.sort_key_val(
      logits + z, iota, dimension=-1, is_stable=False)
  return perm[::-1]

def make_masked_msa(key, batch, opt=None, epsilon=1e-6):
  """Create data for BERT on raw MSA."""
  # Add a random amino acid uniformly.
  
  if opt is None:
    opt = dict(replace_fraction=0.15,
               uniform_prob=0.1,
               profile_prob=0.1, 
               same_prob=0.1)
  
  random_aa = jnp.array([0.05] * 20 + [0., 0.], dtype=jnp.float32)

  msa_one_hot = st_one_hot(batch['msa'], 22)
  categorical_probs = (
      opt["uniform_prob"] * random_aa +
      opt["profile_prob"] * batch['msa_profile'] +
      opt["same_prob"] * msa_one_hot)

  # Put all remaining probability on [MASK] which is a new column.
  pad_shapes = [[0, 0] for _ in range(len(categorical_probs.shape))]
  pad_shapes[-1][1] = 1
  mask_prob = 1. - opt["profile_prob"] - opt["same_prob"] - opt["uniform_prob"]
  categorical_probs = jnp.pad(categorical_probs, pad_shapes, constant_values=mask_prob)
  sh = batch['msa'].shape

  key, mask_subkey, gumbel_subkey = jax.random.split(key, 3)
  uniform = utils.padding_consistent_rng(jax.random.uniform)
  mask_position = uniform(mask_subkey, sh[:2]) < opt["replace_fraction"]
  mask_position *= batch['msa_mask']

  logits = jnp.log(categorical_probs + epsilon)
  bert_msa = gumbel_max_sample(gumbel_subkey, logits)
  msa = jnp.pad(batch['msa'],[[0,0],[0,0],[0,1]])

  bert_msa = jnp.where(mask_position[:,:,None], bert_msa, msa)
  bert_msa *= batch['msa_mask'][:,:,None]

  # Mix real and masked MSA.
  batch['bert_mask'] = mask_position
  batch['true_msa'] = msa
  batch['msa'] = bert_msa

def st_one_hot(x, num_classes, axis=-1):
  '''one_hot with straight-through, to allow gradients to pass'''
  y = jax.nn.one_hot(x.argmax(axis),
    num_classes=num_classes,
    dtype=x.dtype)
  return jax.lax.stop_gradient(y - x) + x

def nearest_neighbor_clusters(batch, gap_agreement_weight=0.):
  """Assign each extra MSA sequence to its nearest neighbor in sampled MSA."""

  # Determine how much weight we assign to each agreement.  In theory, we could
  # use a full blosum matrix here, but right now let's just down-weight gap
  # agreement because it could be spurious.
  # Never put weight on agreeing on BERT mask.

  weights = jnp.array([1.] * 21 + [gap_agreement_weight] + [0.], dtype=jnp.float32)

  msa_mask = batch['msa_mask']
  msa_one_hot = st_one_hot(batch['msa'],23)

  extra_mask = batch['extra_msa_mask']
  extra_one_hot = st_one_hot(batch['extra_msa'],23)

  msa_one_hot_masked = msa_mask[:, :, None] * msa_one_hot
  extra_one_hot_masked = extra_mask[:, :, None] * extra_one_hot

  agreement = jnp.einsum('mrc, nrc->nm', extra_one_hot_masked,
                         weights * msa_one_hot_masked)

  cluster_assignment = jax.nn.softmax(1e3 * agreement, axis=0)
  cluster_assignment *= jnp.einsum('mr, nr->mn', msa_mask, extra_mask)

  cluster_count = jnp.sum(cluster_assignment, axis=-1)
  cluster_count += 1.  # We always include the sequence itself.

  msa_sum = jnp.einsum('nm, mrc->nrc', cluster_assignment, extra_one_hot_masked)
  msa_sum += msa_one_hot_masked

  cluster_profile = msa_sum / cluster_count[:, None, None]

  extra_deletion_matrix = batch['extra_deletion_matrix']
  deletion_matrix = batch['deletion_matrix']

  del_sum = jnp.einsum('nm, mc->nc', cluster_assignment,
                       extra_mask * extra_deletion_matrix)
  del_sum += deletion_matrix  # Original sequence.
  cluster_deletion_mean = del_sum / cluster_count[:, None]

  batch['cluster_profile'] = cluster_profile
  batch['cluster_deletion_mean'] = cluster_deletion_mean

def create_msa_feat(batch):
  """Create and concatenate MSA features."""
  msa = batch['msa']
  deletion_matrix = batch['deletion_matrix']
  
  c_msa = batch.get("cluster_profile",st_one_hot(msa,23))
  c_deletion_matrix = batch.get("cluster_deletion_mean",deletion_matrix)
  
  has_deletion = jnp.clip(deletion_matrix, 0., 1.)[..., None]
  deletion_value = (jnp.arctan(deletion_matrix / 3.) * (2. / jnp.pi))[..., None]
  deletion_mean_value = (jnp.arctan(c_deletion_matrix / 3.) * (2. / jnp.pi))[..., None]
  msa_feat = [
      msa,
      has_deletion,
      deletion_value,
      c_msa,
      deletion_mean_value
  ]
  batch["msa_feat"] = jnp.concatenate(msa_feat, axis=-1)

def create_extra_msa_feature(batch, num_extra_msa):
  """Expand extra_msa into 1hot and concat with other extra msa features.

  We do this as late as possible as the one_hot extra msa can be very large.

  Args:
    batch: a dictionary with the following keys:
     * 'extra_msa': [num_seq, num_res] MSA that wasn't selected as a cluster
       centre. Note - This isn't one-hotted.
     * 'extra_deletion_matrix': [num_seq, num_res] Number of deletions at given
        position.
    num_extra_msa: Number of extra msa to use.

  Returns:
    Concatenated tensor of extra MSA features.
  """
  # 23 = 20 amino acids + 'X' for unknown + gap + bert mask
  extra_msa = batch['extra_msa'][:num_extra_msa]
  deletion_matrix = batch['extra_deletion_matrix'][:num_extra_msa]
  msa_1hot = extra_msa
  has_deletion = jnp.clip(deletion_matrix, 0., 1.)[..., None]
  deletion_value = (jnp.arctan(deletion_matrix / 3.) * (2. / jnp.pi))[..., None]

  batch["extra_msa_feat"] = jnp.concatenate([msa_1hot, has_deletion, deletion_value], axis=-1)
  batch["extra_msa_mask"] = batch['extra_msa_mask'][:num_extra_msa]

def sample_msa(key, batch, max_seq):
  """Sample MSA randomly, remaining sequences are stored as `extra_*`.
  Args:
    key: safe key for random number generation.
    batch: batch to sample msa from.
    max_seq: number of sequences to sample.
  Returns:
    Protein with sampled msa.
  """
  # Sample uniformly among sequences with at least one non-masked position.
  logits = (jnp.clip(jnp.sum(batch['msa_mask'], axis=-1), 0., 1.) - 1.) * 1e6
  # The cluster_bias_mask can be used to preserve the first row (target
  # sequence) for each chain, for example.
  if 'cluster_bias_mask' not in batch:
    cluster_bias_mask = jnp.pad(jnp.zeros(batch['msa'].shape[0] - 1), (1, 0), constant_values=1.)
  else:
    cluster_bias_mask = batch['cluster_bias_mask']

  logits += cluster_bias_mask * 1e6
  index_order = gumbel_argsort_sample_idx(key, logits)
  sel_idx = index_order[:max_seq]
  extra_idx = index_order[max_seq:]

  for k in ['msa', 'deletion_matrix', 'msa_mask']:
    batch['extra_' + k] = batch[k][extra_idx]
    batch[k] = batch[k][sel_idx]

def make_msa_profile(batch):
  """Compute the MSA profile."""
  # Compute the profile for every residue (over all MSA sequences).
  msa_1hot = st_one_hot(batch['msa'],22)
  batch["msa_profile"] = utils.mask_mean(batch['msa_mask'][:,:,None], msa_1hot, axis=0)

def pad_msa_N(batch, num_msa, num_extra_msa):
  if batch["msa"].shape[0] < (num_msa + num_extra_msa):
    N = num_msa + num_extra_msa - batch["msa"].shape[0]
    batch["msa"] = jnp.pad(batch["msa"],[[0,N],[0,0],[0,0]])
    batch["msa_mask"] = jnp.pad(batch["msa_mask"],[[0,N],[0,0]])
    batch["deletion_matrix"] = jnp.pad(batch["deletion_matrix"],[[0,N],[0,0]])

def pad_msa_A(batch, alphabet=23):
  for k in ["msa","extra_msa"]:
    A = alphabet - batch[k].shape[-1]
    batch[k] = jnp.pad(batch[k],[[0,0],[0,0],[0,A]])

def make_msa_feats(batch, key, 
                   num_msa=512, 
                   num_extra_msa=1024,
                   use_mlm=False, mlm_opt=None,
                   use_cluster_profile=False):

  key, sample_key, mask_key = jax.random.split(key,3)

  # mask msa based on seq_mask
  batch["msa_mask"] = jnp.where(batch["seq_mask"],batch["msa_mask"], 0)

  # set num_msa/num_extra_msa
  num = batch["msa"].shape[0]
  if num < num_msa:
    (num_msa,num_extra_msa) = (num,1)
  elif num < (num_msa + num_extra_msa):
    num_extra_msa = num - num_msa

  pad_msa_N(batch, num_msa, num_extra_msa)
  # (msa, msa_mask, deletion_matrix) → (msa, msa_mask, deletion_matrix) 

  make_msa_profile(batch)
  # (msa, msa_mask) → (msa_profile)
  
  sample_msa(sample_key, batch, num_msa)
  # (msa, deletion_matrix, msa_mask)→ ([extra_]msa, [extra_]deletion_matrix, [extra_]msa_mask)
  
  if use_mlm: make_masked_msa(mask_key, batch, mlm_opt)
  # (msa, msa_mask, msa_profile) → (msa, true_msa, bert_mask)
  
  pad_msa_A(batch)
  # ([extra_]msa) → ([extra_]msa)

  if use_cluster_profile: nearest_neighbor_clusters(batch)
  # ([extra_]msa, [extra_]msa_mask, [extra_]deletion_matrix) → (cluster_profile, cluster_deletion_mean)

  create_msa_feat(batch)
  # (msa, deletion_matrix, cluster_profile, cluster_deletion_mean) → (msa_feat)

  create_extra_msa_feature(batch, num_extra_msa)
  # (extra_msa, extra_deletion_matrix, extra_msa_mask) → (extra_msa_feat, extra_msa_mask)

  return batch