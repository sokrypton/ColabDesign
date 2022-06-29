import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.af.alphafold.model import data, config, model

from colabdesign.af.prep import _af_prep
from colabdesign.af.loss import _af_loss
from colabdesign.af.utils import _af_utils
from colabdesign.af.design import _af_design

################################################################
# MK_DESIGN_MODEL - initialize model, and put it all together
################################################################

class mk_afdesign_model(_af_prep, _af_loss, _af_design, _af_utils):
  def __init__(self, protocol="fixbb", num_seq=1,
               num_models=1, model_sample=True,
               recycle_mode="average", num_recycles=0,
               use_templates=False, data_dir=".",
               debug=False, loss_callback=None):
    
    # decide if templates should be used
    if protocol == "binder": use_templates = True

    self.protocol = protocol

    self._loss_callback = loss_callback
    self._num = num_seq    
    self._args = {"use_templates":use_templates,
                  "recycle_mode": recycle_mode,
                  "model_sample": model_sample,
                  "debug":debug, "repeat": False}
    
    self._opt = {"dropout":True, "lr":1.0,
                 "recycles":num_recycles, "models":num_models,
                 "temp":1.0, "soft":1.0, "hard":1.0, "bias":0.0,
                 "con":      {"num":2, "cutoff":14.0, "binary":False, "seqsep":9},
                 "i_con":    {"num":1, "cutoff":20.0, "binary":False},                 
                 "template": {"aatype":21, "dropout":0.15},
                 "weights":  {"helix":0.0, "plddt":0.01, "pae":0.01}}

    #############################
    # CONFIG
    #############################
    # decide which condig to use configs to use
    if use_templates:
      model_name = "model_1_ptm"
      self._opt["models"] = min(num_models, 2)
    else:
      model_name = "model_3_ptm"
    
    cfg = config.model_config(model_name)
    cfg.model.global_config.use_remat = True  
    cfg.model.global_config.subbatch_size = None    
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

    ##################################
    # MODEL
    ##################################
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
    self._grad_fn, self._fn = [jax.jit(x) for x in self._get_model()]

    #####################################
    # set protocol specific functions
    #####################################
    protocols = ["fixbb","hallucination","binder","partial"]
    prep_fns = [self._prep_fixbb, self._prep_hallucination, self._prep_binder, self._prep_partial]
    loss_fns = [self._loss_fixbb, self._loss_hallucination, self._loss_binder, self._loss_partial]
    self.prep_inputs = dict(zip(protocols, prep_fns))[self.protocol]
    self._get_loss = dict(zip(protocols, loss_fns))[self.protocol]