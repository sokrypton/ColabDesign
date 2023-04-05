import os,jax
# disable triton_gemm for jax versions > 0.3
if int(jax.__version__.split(".")[1]) > 3:
  os.environ["XLA_FLAGS"] = "--xla_gpu_enable_triton_gemm=false"

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from colabdesign.shared.utils import clear_mem
from .model import mk_mpnn_model