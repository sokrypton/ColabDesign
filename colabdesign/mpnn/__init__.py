import os
os.environ["XLA_FLAGS"] = "--xla_gpu_enable_triton_gemm=false"

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from colabdesign.shared.utils import clear_mem
from .model import mk_mpnn_model