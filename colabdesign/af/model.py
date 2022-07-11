import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.af.alphafold.model import data, config, model

from colabdesign.shared.model import design_model
from colabdesign.shared.utils import Key

from colabdesign.af.prep   import _af_prep
from colabdesign.af.loss   import _af_loss, get_plddt, get_pae
from colabdesign.af.utils  import _af_utils
from colabdesign.af.design import _af_design
from colabdesign.af.inputs import _af_inputs, update_seq, update_aatype 

################################################################
# MK_DESIGN_MODEL - initialize model, and put it all together
################################################################

class mk_af_model(design_model, _af_inputs, _af_loss, _af_prep, _af_design, _af_utils):
  def __init__(self, protocol="fixbb", num_seq=1,
               num_models=1, sample_models=True,
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
                  "debug":debug, "repeat": False}
    
    self.opt = {"dropout":True, "lr":1.0,
                "recycles":num_recycles, "models":num_models, "sample_models":sample_models,
                "temp":1.0, "soft":0.0, "hard":0.0, "bias":0.0, "alpha":2.0,
                "con":      {"num":2, "cutoff":14.0, "binary":False, "seqsep":9},
                "i_con":    {"num":1, "cutoff":20.0, "binary":False},                 
                "template": {"aatype":21, "dropout":0.15},
                "weights":  {"helix":0.0, "plddt":0.01, "pae":0.01}}
    
    self.params = {}

    #############################
    # configure AlphaFold
    #############################
    # decide which condig to use configs to use
    if use_templates:
      model_name = "model_1_ptm"
      self.opt["models"] = min(num_models, 2)
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
      cfg.data.eval.max_templates = 0
      cfg.data.eval.max_msa_clusters = num_seq
    cfg.data.common.max_extra_msa = 1
    cfg.data.eval.masked_msa_replace_fraction = 0

    # number of recycles
    assert recycle_mode in ["average","add_prev","backprop","last","sample"]
    if recycle_mode == "average": num_recycles = 0
    cfg.data.common.num_recycle = 0      # for feature processing
    cfg.model.num_recycle = num_recycles # for model configuration
    self._config = cfg
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
    idx = ["fixbb","hallucination","binder","partial"].index(self.protocol)
    self.prep_inputs = [self._prep_fixbb, self._prep_hallucination, self._prep_binder, self._prep_partial][idx]
    self._get_loss   = [self._loss_fixbb, self._loss_hallucination, self._loss_binder, self._loss_partial][idx]

  def _get_model(self, callback=None):

    # setup function to get gradients
    def _model(params, model_params, inputs, key, opt):

      aux = {}
      key = Key(key=key).get

      #######################################################################
      # INPUTS
      #######################################################################

      # get sequence
      seq = self._get_seq(params, opt, aux, key())
            
      # update sequence features
      update_seq(seq["pseudo"], inputs)
      
      # update amino acid sidechain identity
      B,L = inputs["aatype"].shape[:2]
      aatype = jax.nn.one_hot(seq["pseudo"][0].argmax(-1),21)
      update_aatype(jnp.broadcast_to(aatype,(B,L,21)), inputs)
      
      # update template features
      if self._args["use_templates"]:
        self._update_template(inputs, opt, key())
      
      # set dropout
      inputs["dropout_scale"] = jnp.array([opt["dropout"]]).astype(float)
      
      # decide number of recycles to do
      if self._args["recycle_mode"] in ["last","sample"]:
        inputs["num_iter_recycling"] = jnp.array([opt["recycles"]])

      #######################################################################
      # OUTPUTS
      #######################################################################
      outputs = self._runner.apply(model_params, key(), inputs)

      # add aux outputs
      aux.update({"final_atom_positions":outputs["structure_module"]["final_atom_positions"],
                  "final_atom_mask":outputs["structure_module"]["final_atom_mask"],
                  "plddt":get_plddt(outputs),"pae":get_pae(outputs),
                  "prev":outputs["prev"]})

      if self._args["debug"]:
        aux["debug"] = {"inputs":inputs, "outputs":outputs, "opt":opt}

      #######################################################################
      # LOSS
      #######################################################################
      aux["losses"] = {}
      self._get_loss(inputs=inputs, outputs=outputs, opt=opt, aux=aux)

      if self._loss_callback is not None:
        aux["losses"].update(self._loss_callback(inputs, outputs, opt))

      # weighted loss
      w = opt["weights"]
      loss = sum([v*w[k] if k in w else v for k,v in aux["losses"].items()])
            
      return loss, aux
    
    return jax.value_and_grad(_model, has_aux=True, argnums=0), _model