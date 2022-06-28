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

"""Code for constructing the model."""
from typing import Any, Mapping, Optional, Union

from absl import logging
from colabdesign.af.alphafold.common import confidence
from colabdesign.af.alphafold.model import features
from colabdesign.af.alphafold.model import modules
import haiku as hk
import jax
import ml_collections
import numpy as np
import tensorflow.compat.v1 as tf
import tree


def get_confidence_metrics(
    prediction_result: Mapping[str, Any]) -> Mapping[str, Any]:
  """Post processes prediction_result to get confidence metrics."""

  confidence_metrics = {}
  confidence_metrics['plddt'] = confidence.compute_plddt(
      prediction_result['predicted_lddt']['logits'])
  if 'predicted_aligned_error' in prediction_result:
    confidence_metrics.update(confidence.compute_predicted_aligned_error(
        prediction_result['predicted_aligned_error']['logits'],
        prediction_result['predicted_aligned_error']['breaks']))
    confidence_metrics['ptm'] = confidence.predicted_tm_score(
        prediction_result['predicted_aligned_error']['logits'],
        prediction_result['predicted_aligned_error']['breaks'])

  return confidence_metrics


class RunModel:
  """Container for JAX model."""

  def __init__(self,
               config: ml_collections.ConfigDict,
               params: Optional[Mapping[str, Mapping[str, np.ndarray]]] = None,
               is_training=True,
               return_representations=True,
               recycle_mode=None):
    self.config = config
    self.params = params

    
    self.mode = recycle_mode
    if self.mode is None:
      self.mode = []
      # backward compatibility
      if self.config.model.add_prev:
        self.mode.append("add_prev")
      if self.config.model.backprop_recycle:
        self.mode.append("backprop")

    def _forward_fn(batch):
      model = modules.AlphaFold(self.config.model)
      return model(
          batch,
          is_training=is_training,
          compute_loss=False,
          ensemble_representations=False,
          return_representations=return_representations)
    
    self.init = jax.jit(hk.transform(_forward_fn).init)
    self.apply_fn = jax.jit(hk.transform(_forward_fn).apply)
    
    def apply(params, key, feat):
      
      if "prev" in feat:
        prev = feat["prev"]      
      else:
        L = feat['aatype'].shape[1]
        prev = {'prev_msa_first_row': np.zeros([L,256]),
                'prev_pair': np.zeros([L,L,128])}
        if self.config.model.use_struct:
          prev['prev_pos'] = np.zeros([L,37,3])
        else:
          prev['prev_dgram'] = np.zeros([L,L,64])

      ################################
      # decide how to run recycles
      ################################
      if "num_iter_recycling" in feat:
        # use while_loop()
        num_recycles = feat.pop("num_iter_recycling")[0]
        def body(x):
          i,prev,key = x
          key, sub_key = jax.random.split(key)
          feat["prev"] = prev
          prev = self.apply_fn(params, sub_key, feat)["prev"]
          prev = jax.lax.stop_gradient(prev)
          return (i+1, prev, key)

        init = (0,prev,key)
        _, feat["prev"], key = jax.lax.while_loop(lambda x: x[0] < num_recycles, body, init)        
        key, sub_key = jax.random.split(key)
        results = self.apply_fn(params, sub_key, feat)

      elif self.config.model.num_recycle:
        # use scan()
        def loop(prev, sub_key):
          feat["prev"] = prev
          results = self.apply_fn(params, sub_key, feat)
          prev = results["prev"]
          if "backprop" not in self.mode:
              prev = jax.lax.stop_gradient(prev)
          return prev, results

        keys = jax.random.split(key, self.config.model.num_recycle + 1)
        _, o = jax.lax.scan(loop, prev, keys)
        results = jax.tree_map(lambda x:x[-1], o)

        if "add_prev" in self.mode:
          for k in ["distogram","predicted_lddt","predicted_aligned_error"]:
            if k in results:
              results[k]["logits"] = o[k]["logits"].mean(0)
      
      else:
        # single pass
        results = self.apply_fn(params, key, feat)
      
      return results
    
    self.apply = jax.jit(apply)

  def init_params(self, feat: features.FeatureDict, random_seed: int = 0):
    """Initializes the model parameters.

    If none were provided when this class was instantiated then the parameters
    are randomly initialized.

    Args:
      feat: A dictionary of NumPy feature arrays as output by
        RunModel.process_features.
      random_seed: A random seed to use to initialize the parameters if none
        were set when this class was initialized.
    """
    if not self.params:
      # Init params randomly.
      rng = jax.random.PRNGKey(random_seed)
      self.params = hk.data_structures.to_mutable_dict(
          self.init(rng, feat))
      logging.warning('Initialized parameters randomly')

  def process_features(
      self,
      raw_features: Union[tf.train.Example, features.FeatureDict],
      random_seed: int) -> features.FeatureDict:
    """Processes features to prepare for feeding them into the model.

    Args:
      raw_features: The output of the data pipeline either as a dict of NumPy
        arrays or as a tf.train.Example.
      random_seed: The random seed to use when processing the features.

    Returns:
      A dict of NumPy feature arrays suitable for feeding into the model.
    """
    if isinstance(raw_features, dict):
      return features.np_example_to_features(
          np_example=raw_features,
          config=self.config,
          random_seed=random_seed)
    else:
      return features.tf_example_to_features(
          tf_example=raw_features,
          config=self.config,
          random_seed=random_seed)

  def eval_shape(self, feat: features.FeatureDict) -> jax.ShapeDtypeStruct:
    self.init_params(feat)
    logging.info('Running eval_shape with shape(feat) = %s', tree.map_structure(lambda x: x.shape, feat))
    shape = jax.eval_shape(self.apply, self.params, jax.random.PRNGKey(0), feat)
    logging.info('Output shape was %s', shape)
    return shape

  def predict(self, feat: features.FeatureDict) -> Mapping[str, Any]:
    """Makes a prediction by inferencing the model on the provided features.

    Args:
      feat: A dictionary of NumPy feature arrays as output by
        RunModel.process_features.

    Returns:
      A dictionary of model outputs.
    """
    self.init_params(feat)
    logging.info('Running predict with shape(feat) = %s', tree.map_structure(lambda x: x.shape, feat))

    result = self.apply(self.params, jax.random.PRNGKey(0), feat)
    if self.config.model.use_struct:
      result.update(get_confidence_metrics(result))
        
    logging.info('Output shape was %s', tree.map_structure(lambda x: x.shape, result))
    return result
