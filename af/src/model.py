import jax
import jax.numpy as jnp
import numpy as np

from alphafold.model import data, config, model

from af.src.prep import _af_prep
from af.src.loss import _af_loss
from af.src.utils import _af_utils
from af.src.design import _af_design

################################################################
# MK_DESIGN_MODEL - initialize model, and put it all together
################################################################

class mk_afdesign_model(_af_prep, _af_loss, _af_design, _af_utils):
  def __init__(self, protocol="fixbb", num_seq=1,
               num_models=1, model_sample=True,
               recycle_mode="average", num_recycles=0,
               use_templates=False, use_pssm=False, data_dir=".",
               debug=False, use_struct=True):
    
    # decide if templates should be used
    if protocol == "binder": use_templates = True

    self.protocol = protocol
    self.use_struct = use_struct
    self.args = {"num_seq":num_seq, "use_templates":use_templates,
                 "recycle_mode": recycle_mode,
                 "model_sample":model_sample,
                 "use_pssm":use_pssm, "debug":debug,
                 "repeat": False}
    
    self._default_opt = {
                         "temp":1.0, "soft":1.0, "hard":1.0, "gumbel":False,
                         "dropout":True, "dropout_scale":1.0,
                         "recycles":num_recycles, "models":num_models,
                         "con":  {"num":2, "cutoff":14.0, "seqsep":9, "binary":False, "entropy":True},
                         "i_con":{"num":1, "cutoff":20.0,             "binary":False, "entropy":True},
                         "bias":0.0, "template_aatype":21, "template_dropout":0.15,
                         # default weights
                         "weights:{"msa_ent":0.0, "helix":0.0, "plddt":0.01, "pae":0.01}
                        }

    # setup which model configs to use
    if use_templates:
      model_name = "model_1_ptm"
      self._default_opt["models"] = min(num_models, 2)
    else:
      model_name = "model_3_ptm"
    
    cfg = config.model_config(model_name)

    # enable checkpointing
    cfg.model.global_config.use_remat = True  

    # subbatch_size / chunking
    cfg.model.global_config.subbatch_size = None
    
    # enable/disable structure module
    cfg.model.use_struct = use_struct

    # number of sequences
    if use_templates:
      cfg.data.eval.max_templates = 1
      cfg.data.eval.max_msa_clusters = num_seq + 1
    else:
      cfg.data.eval.max_msa_clusters = num_seq

    cfg.data.common.max_extra_msa = 1
    cfg.data.eval.masked_msa_replace_fraction = 0

    # number of recycles (recycles are now managed in AfDesign)
    assert recycle_mode in ["average","add_prev","backprop","last","sample"]
    if recycle_mode == "average": num_recycles = 0
    
    cfg.data.common.num_recycle = 0      # for feature processing
    cfg.model.num_recycle = num_recycles # for model configuration

    self._config = cfg

    # setup model
    self._model_params = [data.get_model_haiku_params(model_name=model_name, data_dir=data_dir)]
    self._runner = model.RunModel(self._config, self._model_params[0],
                                  is_training=True, recycle_mode=recycle_mode)

    # load the other model_params
    if use_templates:
      model_names = ["model_2_ptm"]
    else:
      model_names = ["model_1_ptm","model_2_ptm","model_4_ptm","model_5_ptm"]
    for model_name in model_names:
      params = data.get_model_haiku_params(model_name=model_name, data_dir=data_dir)
      self._model_params.append({k: params[k] for k in self._runner.params.keys()})

    # define gradient function
    self._grad_fn, self._fn = [jax.jit(x) for x in self._get_fn()]

    # define input function
    if protocol == "fixbb":           self.prep_inputs = self._prep_fixbb
    if protocol == "hallucination":   self.prep_inputs = self._prep_hallucination
    if protocol == "binder":          self.prep_inputs = self._prep_binder
    if protocol == "partial":         self.prep_inputs = self._prep_partial
