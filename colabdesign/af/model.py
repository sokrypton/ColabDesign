import os
import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.af.alphafold.model import data, config, model

from colabdesign.shared.model import design_model
from colabdesign.shared.utils import Key

from colabdesign.af.prep   import _af_prep
from colabdesign.af.loss   import _af_loss, get_plddt, get_pae, get_contact_map
from colabdesign.af.utils  import _af_utils
from colabdesign.af.design import _af_design
from colabdesign.af.inputs import _af_inputs, update_seq, update_aatype, crop_feat

################################################################
# MK_DESIGN_MODEL - initialize model, and put it all together
################################################################

class mk_af_model(design_model, _af_inputs, _af_loss, _af_prep, _af_design, _af_utils):
  def __init__(self, protocol="fixbb", num_seq=1,
               num_models=1, sample_models=True,
               recycle_mode="average", num_recycles=0,
               use_templates=False, best_metric="loss",
               crop_len=None, crop_mode="pair",
               debug=False, use_alphafold=True, use_openfold=False,
               loss_callback=None, data_dir="."):
    
    assert protocol in ["fixbb","hallucination","binder","partial"]
    assert recycle_mode in ["average","add_prev","backprop","last","sample"]
    assert crop_mode in ["slide","roll","pair"]
    
    # decide if templates should be used
    if protocol == "binder": use_templates = True

    self.protocol = protocol
    self._loss_callback = loss_callback
    self._num = num_seq
    self._args = {"use_templates":use_templates,
                  "recycle_mode":recycle_mode,
                  "debug":debug,
                  "repeat":False, "homooligomer":False, "copies":1,
                  "best_metric":best_metric,
                  'use_alphafold':use_alphafold, 'use_openfold':use_openfold,
                  "crop_len":crop_len,"crop_mode":crop_mode}
    
    self.opt = {"dropout":True, "lr":1.0, "use_pssm":False,
                "recycles":num_recycles, "models":num_models, "sample_models":sample_models,
                "temp":1.0, "soft":0.0, "hard":0.0, "bias":0.0, "alpha":2.0,
                "con":      {"num":2, "cutoff":14.0, "binary":False, "seqsep":9},
                "i_con":    {"num":1, "cutoff":20.0, "binary":False},                 
                "template": {"aatype":21, "dropout":0.0},
                "weights":  {"helix":0.0, "plddt":0.01, "pae":0.01},
                "cmap_cutoff":10.0}
    
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
    if recycle_mode == "average": num_recycles = 0
    cfg.data.common.num_recycle = 0      # for feature processing
    cfg.model.num_recycle = num_recycles # for model configuration

    # setup model
    self._cfg = cfg 

    # load model_params
    model_names = []
    if use_templates:
      model_names += [f"model_{k}_ptm" for k in [1,2]]
      model_names += [f"openfold_model_ptm_{k}" for k in [1,2]]    
    else:
      model_names += [f"model_{k}_ptm" for k in [1,2,3,4,5]]
      model_names += [f"openfold_model_ptm_{k}" for k in [1,2]] + ["openfold_model_no_templ_ptm_1"]

    self._model_params, self._model_names = [],[]
    for model_name in model_names:
      params = data.get_model_haiku_params(model_name=model_name, data_dir=data_dir)
      if params is not None:
        if not use_templates:
          params = {k:v for k,v in params.items() if "template" not in k}
        self._model_params.append(params)
        self._model_names.append(model_name)

    #####################################
    # set protocol specific functions
    #####################################
    idx = ["fixbb","hallucination","binder","partial"].index(self.protocol)
    self.prep_inputs = [self._prep_fixbb, self._prep_hallucination, self._prep_binder, self._prep_partial][idx]
    self._get_loss   = [self._loss_fixbb, self._loss_hallucination, self._loss_binder, self._loss_partial][idx]

  def _get_model(self, cfg, callback=None):

    runner = model.RunModel(cfg, is_training=True, recycle_mode=self._args["recycle_mode"])

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
      pssm = jnp.where(opt["use_pssm"], seq["pssm"], seq["pseudo"])
      update_seq(seq["pseudo"], inputs, seq_pssm=pssm)
      
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

      # batch
      batch = self._batch if hasattr(self,"_batch") else None

      # crop inputs
      if opt["crop_pos"].shape[0] < L:
        inputs = crop_feat(inputs, opt["crop_pos"], self._cfg, add_batch=True)    
        batch = crop_feat(batch, opt["crop_pos"], self._cfg, add_batch=False)

      #######################################################################
      # OUTPUTS
      #######################################################################
      outputs = runner.apply(model_params, key(), inputs)

      # add aux outputs
      aux.update({"atom_positions":outputs["structure_module"]["final_atom_positions"],
                  "atom_mask":outputs["structure_module"]["final_atom_mask"],                  
                  "residue_index":inputs["residue_index"][0], "aatype":inputs["aatype"][0],
                  "plddt":get_plddt(outputs),"pae":get_pae(outputs),
                  "cmap":get_contact_map(outputs, opt["cmap_cutoff"])})

      # experimental
      # crop outputs (TODO)
      if opt["crop_pos"].shape[0] < L:
        p = opt["crop_pos"]
        aux["cmap"] = jnp.zeros((L,L)).at[p[:,None],p[None,:]].set(aux["cmap"])
        aux["pae"] = jnp.full((L,L),jnp.nan).at[p[:,None],p[None,:]].set(aux["pae"])

      if self._args["recycle_mode"] == "average": aux["prev"] = outputs["prev"]
      if self._args["debug"]: aux["debug"] = {"inputs":inputs, "outputs":outputs, "opt":opt}

      #######################################################################
      # LOSS
      #######################################################################
      aux["losses"] = {}
      self._get_loss(inputs=inputs, outputs=outputs, opt=opt, aux=aux, batch=batch)

      if self._loss_callback is not None:
        aux["losses"].update(self._loss_callback(inputs, outputs, opt))

      # weighted loss
      w = opt["weights"]
      loss = sum([v * w[k] if k in w else v for k,v in aux["losses"].items()])
            
      return loss, aux
    
    return {"grad_fn":jax.jit(jax.value_and_grad(_model, has_aux=True, argnums=0)),
            "fn":jax.jit(_model)}