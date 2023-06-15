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
               use_multimer=False):

    self.config = config    
    def _forward_fn(batch):
      if use_multimer:
        return modules_multimer.AlphaFold(self.config.model)(batch)
      else:
        return modules.AlphaFold(self.config.model)(batch)    
        
    self.apply = jax.jit(hk.transform(_forward_fn).apply)