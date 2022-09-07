import os
import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.af.alphafold.model import data, config, model, all_atom

from colabdesign.shared.model import design_model
from colabdesign.shared.utils import Key

from colabdesign.af.prep   import _af_prep
from colabdesign.af.loss   import _af_loss, get_plddt, get_pae, get_contact_map, get_ptm, get_seq_ent_loss, get_mlm_loss
from colabdesign.af.utils  import _af_utils
from colabdesign.af.design import _af_design
from colabdesign.af.inputs import _af_inputs, update_seq, update_aatype
from colabdesign.af.crop   import _af_crop, crop_feat

################################################################
# MK_DESIGN_MODEL - initialize model, and put it all together
################################################################

class mk_af_model(design_model, _af_inputs, _af_loss, _af_prep, _af_design, _af_utils, _af_crop):
  def __init__(self, protocol="fixbb", num_seq=1,
               num_models=1, sample_models=True,
               recycle_mode="last", num_recycles=0,
               use_templates=False, best_metric="loss",
               model_names=None,
               use_openfold=False, use_alphafold=True,
               use_multimer=False,
               use_mlm=False,
               use_crop=False, crop_len=None, crop_mode="slide",
               debug=False, loss_callback=None, data_dir="."):
    
    assert protocol in ["fixbb","hallucination","binder","partial"]
    assert recycle_mode in ["average","first","last","sample","add_prev","backprop"]
    assert crop_mode in ["slide","roll","pair","dist"]
    
    # decide if templates should be used
    if protocol == "binder": use_templates = True

    self.protocol = protocol
    self._loss_callback = loss_callback
    self._num = num_seq
    self._args = {"use_templates":use_templates, "use_multimer":use_multimer,
                  "recycle_mode":recycle_mode, "use_mlm": use_mlm,
                  "debug":debug,
                  "repeat":False, "homooligomer":False, "copies":1,
                  "best_metric":best_metric,
                  "use_crop":use_crop, "crop_len":crop_len, "crop_mode":crop_mode,
                  "models":None}

    self.opt = {"dropout":True, "lr":1.0, "use_pssm":False, "mlm_dropout":0.05,
                "num_recycles":num_recycles, "num_models":num_models, "sample_models":sample_models,
                "temp":1.0, "soft":0.0, "hard":0.0, "bias":0.0, "alpha":2.0,
                "con":      {"num":2, "cutoff":14.0, "binary":False, "seqsep":9, "num_pos":float("inf")},
                "i_con":    {"num":1, "cutoff":21.6875, "binary":False, "num_pos":float("inf")},
                "template": {"dropout":0.0, "rm_ic":False, "rm_seq":True, "rm_sc":True},                
                "weights":  {"seq_ent":0.0, "plddt":0.0, "pae":0.0, "exp_res":0.0},
                "cmap_cutoff": 10.0, "fape_cutoff":10.0}

    if self._args["use_mlm"]:
      self.opt["weights"]["mlm"] = 0.1
    
    self._params = {}
    self._inputs = {}

    #############################
    # configure AlphaFold
    #############################
    if use_multimer:
      cfg = config.model_config("model_1_multimer")
    else:
      cfg = config.model_config("model_1_ptm" if use_templates else "model_3_ptm")
    if recycle_mode in ["average","first","last","sample"]: num_recycles = 0
    cfg.model.num_recycle = num_recycles
    cfg.model.global_config.use_remat = True

    # setup model
    self._cfg = cfg 

    # load model_params
    if model_names is None:
      model_names = []
      if use_multimer:
        model_names += [f"model_{k}_multimer_v2" for k in [1,2,3,4,5]]
      else:
        if use_templates:
          if use_alphafold: model_names += [f"model_{k}_ptm" for k in [1,2]]
          if use_openfold:  model_names += [f"openfold_model_ptm_{k}" for k in [1,2]]    
        else:
          if use_alphafold: model_names += [f"model_{k}_ptm" for k in [1,2,3,4,5]]
          if use_openfold:  model_names += [f"openfold_model_ptm_{k}" for k in [1,2]] + ["openfold_model_no_templ_ptm_1"]

    self._model_params, self._model_names = [],[]
    for model_name in model_names:
      params = data.get_model_haiku_params(model_name=model_name, data_dir=data_dir)
      if params is not None:
        if not use_multimer and not use_templates:
          params = {k:v for k,v in params.items() if "template" not in k}
        self._model_params.append(params)
        self._model_names.append(model_name)
      else:
        print(f"WARNING: '{model_name}' not found")

    #####################################
    # set protocol specific functions
    #####################################
    idx = ["fixbb","hallucination","binder","partial"].index(self.protocol)
    self.prep_inputs = [self._prep_fixbb, self._prep_hallucination, self._prep_binder, self._prep_partial][idx]
    self._get_loss   = [self._loss_fixbb, self._loss_hallucination, self._loss_binder, self._loss_partial][idx]

  def _get_model(self, cfg, callback=None):

    a = self._args
    runner = model.RunModel(cfg, is_training=True,
                            recycle_mode=a["recycle_mode"],
                            use_multimer=a["use_multimer"])

    # setup function to get gradients
    def _model(params, model_params, inputs, key, opt):

      aux = {}
      key = Key(key=key).get

      #######################################################################
      # INPUTS
      #######################################################################

      L = inputs["aatype"].shape[0]

      # get sequence
      seq = self._get_seq(inputs, params, opt, aux, key())
            
      # update sequence features
      pssm = jnp.where(opt["use_pssm"], seq["pssm"], seq["pseudo"])
      if a["use_mlm"]:
        mlm = jax.random.bernoulli(key(), opt["mlm_dropout"], (L,))
        update_seq(seq["pseudo"], inputs, seq_pssm=pssm, mlm=mlm)
      else:
        update_seq(seq["pseudo"], inputs, seq_pssm=pssm)
      
      # update amino acid sidechain identity
      aatype = jax.nn.one_hot(seq["pseudo"][0].argmax(-1),21)
      update_aatype(jnp.broadcast_to(aatype,(L,21)), inputs)
      
      # update template features
      if a["use_templates"]:
        self._update_template(inputs, opt, key())
      inputs["mask_template_interchain"] = opt["template"]["rm_ic"]
      
      # set dropout
      inputs["dropout_scale"] = jnp.array(opt["dropout"], dtype=float)
      
      # experimental - crop inputs
      if a["use_crop"] and opt["crop_pos"].shape[0] < L:
          inputs = crop_feat(inputs, opt["crop_pos"])    

      if "batch" not in inputs:
        inputs["batch"] = None
      
      #######################################################################
      # OUTPUTS
      #######################################################################

      outputs = runner.apply(model_params, key(), inputs)

      # add aux outputs
      aux.update({"atom_positions":outputs["structure_module"]["final_atom_positions"],
                  "atom_mask":outputs["structure_module"]["final_atom_mask"],                  
                  "residue_index":inputs["residue_index"],
                  "aatype":inputs["aatype"],
                  "plddt": get_plddt(outputs),
                  "pae":   get_pae(outputs), 
                  "ptm":   get_ptm(inputs, outputs),
                  "i_ptm": get_ptm(inputs, outputs, interface=True), 
                  "cmap":  get_contact_map(outputs, opt["cmap_cutoff"]),
                  "prev": outputs["prev"]})

      # experimental - uncrop outputs
      if a["use_crop"] and opt["crop_pos"].shape[0] < L:
        p = opt["crop_pos"]
        aux["cmap"] = jnp.zeros((L,L)).at[p[:,None],p[None,:]].set(aux["cmap"])
        aux["pae"] = jnp.full((L,L),jnp.nan).at[p[:,None],p[None,:]].set(aux["pae"])
            
      #######################################################################
      # LOSS
      #######################################################################
      aux["losses"] = {}

      # add protocol specific losses
      self._get_loss(inputs=inputs, outputs=outputs, opt=opt, aux=aux)

      # add user defined losses
      inputs["seq"] = aux["seq"]      
      if self._loss_callback is not None:
        loss_fns = self._loss_callback if isinstance(self._loss_callback,list) else [self._loss_callback]
        for loss_fn in loss_fns:
          aux["losses"].update(loss_fn(inputs, outputs, opt))

      if a["debug"]:
        aux["debug"] = {"inputs":inputs, "outputs":outputs, "opt":opt}

      # sequence entropy loss
      aux["losses"].update(get_seq_ent_loss(inputs, outputs, opt))
      
      # experimental masked-language-modeling
      if a["use_mlm"]:
        aux["losses"].update(get_mlm_loss(outputs, mask=mlm, truth=seq["pssm"]))
  
      # weighted loss
      w = opt["weights"]
      loss = sum([v * w[k] if k in w else v for k,v in aux["losses"].items()])
            
      return loss, aux
    
    return {"grad_fn":jax.jit(jax.value_and_grad(_model, has_aux=True, argnums=0)),
            "fn":jax.jit(_model), "runner":runner}
