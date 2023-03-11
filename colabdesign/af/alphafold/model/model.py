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
from colabdesign.af.alphafold.model import modules
from colabdesign.af.alphafold.model import modules_multimer

import haiku as hk
import jax
import ml_collections
import numpy as np
import tree

class RunModel:
  """Container for JAX model."""

  def __init__(self,
               config: ml_collections.ConfigDict,
               params: Optional[Mapping[str, Mapping[str, np.ndarray]]] = None,
               return_representations=True,
               recycle_mode=None,
               use_multimer=False):

    self.config = config
    self.params = params
    
    self.mode = recycle_mode
    if self.mode is None: self.mode = []

    def _forward_fn(batch):
      if use_multimer:
        model = modules_multimer.AlphaFold(self.config.model)
      else:
        model = modules.AlphaFold(self.config.model)
      return model(
          batch,
          return_representations=return_representations)
    
    self.init = jax.jit(hk.transform(_forward_fn).init)
    self.apply_fn = jax.jit(hk.transform(_forward_fn).apply)
    
    def apply(params, key, feat):
      
      if "prev" in feat:
        prev = feat["prev"]      
      else:
        L = feat['aatype'].shape[0]
        prev = {'prev_msa_first_row': np.zeros([L,256]),
                'prev_pair': np.zeros([L,L,128]),
                'prev_pos': np.zeros([L,37,3])}
        if self.config.global_config.use_dgram:
          prev['prev_dgram'] = np.zeros([L,L,64])
        feat["prev"] = prev

      ################################
      # decide how to run recycles
      ################################
      if self.config.model.num_recycle:
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