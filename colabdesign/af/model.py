import os
import jax
import jax.numpy as jnp
import numpy as np
from inspect import signature

from colabdesign.af.alphafold.model import data, config, model, all_atom

from colabdesign.shared.model import design_model
from colabdesign.shared.utils import Key

from colabdesign.af.prep   import _af_prep
from colabdesign.af.loss   import _af_loss, get_plddt, get_pae, get_ptm
from colabdesign.af.loss   import get_contact_map, get_seq_ent_loss, get_mlm_loss
from colabdesign.af.utils  import _af_utils
from colabdesign.af.design import _af_design
from colabdesign.af.inputs import _af_inputs, update_seq, update_aatype

################################################################
# MK_DESIGN_MODEL - initialize model, and put it all together
################################################################

class mk_af_model(design_model, _af_inputs, _af_loss, _af_prep, _af_design, _af_utils):
  def __init__(self, protocol="fixbb", num_seq=1,
               num_models=1, sample_models=True,
               recycle_mode="last", num_recycles=0,
               use_templates=False, best_metric="loss",
               model_names=None, optimizer="sgd", learning_rate=0.1,
               use_bfloat16=True, use_multimer=False,
               pre_callback=None, post_callback=None,
               pre_design_callback=None, post_design_callback=None,
               loss_callback=None, traj_iter=1, traj_max=10000, debug=False, data_dir="."):
    
    assert protocol in ["fixbb","hallucination","binder","partial"]
    assert recycle_mode in ["average","first","last","sample","add_prev","backprop"]

    # decide if templates should be used
    if protocol == "binder": use_templates = True

    self.protocol = protocol
    self._num = num_seq
    self._args = {"use_templates":use_templates, "use_multimer":use_multimer,
                  "recycle_mode":recycle_mode, "use_mlm": False, "realign": True,
                  "debug":debug, "repeat":False, "homooligomer":False, "copies":1,
                  "optimizer":optimizer, "best_metric":best_metric, 
                  "traj_iter":traj_iter, "traj_max":traj_max,
                  "clear_prev": True, "use_dgram":False, "shuffle_msa":True}

    self.opt = {"dropout":True, "pssm_hard":False, "learning_rate":learning_rate, "norm_seq_grad":True,
                "num_recycles":num_recycles, "num_models":num_models, "sample_models":sample_models,                
                "temp":1.0, "soft":0.0, "hard":0.0, "alpha":2.0,
                "con":      {"num":2, "cutoff":14.0, "binary":False, "seqsep":9, "num_pos":float("inf")},
                "i_con":    {"num":1, "cutoff":21.6875, "binary":False, "num_pos":float("inf")},
                "template": {"dropout":0.0, "rm_ic":False},                
                "weights":  {"seq_ent":0.0, "plddt":0.0, "pae":0.0, "exp_res":0.0, "helix":0.0},
                "fape_cutoff":10.0}

    if self._args["use_mlm"]:
      self.opt["mlm_dropout"] = 0.0
      self.opt["weights"]["mlm"] = 0.1

    self._params = {}
    self._inputs = {}
    self._tmp = {"traj":{"seq":[],"xyz":[],"plddt":[],"pae":[]},
                 "log":[],"best":{}}

    # collect callbacks
    self._callbacks = {"model":{"pre":pre_callback,"post":post_callback,"loss":loss_callback},
                       "design":{"pre":pre_design_callback,"post":post_design_callback}}
    
    for m,n in self._callbacks.items():
      for k,v in n.items():
        if v is None: v = []
        if not isinstance(v,list): v = [v]
        self._callbacks[m][k] = v

    #############################
    # configure AlphaFold
    #############################
    if use_multimer:
      cfg = config.model_config("model_1_multimer")
      self.opt["pssm_hard"] = True
    else:
      cfg = config.model_config("model_1_ptm" if use_templates else "model_3_ptm")
    if recycle_mode in ["average","first","last","sample"]: num_recycles = 0
    cfg.model.num_recycle = num_recycles
    cfg.model.global_config.use_remat = True
    cfg.model.global_config.use_dgram = self._args["use_dgram"]
    cfg.model.global_config.bfloat16 = use_bfloat16

    # setup model
    self._cfg = cfg 

    # load model_params
    if model_names is None:
      model_names = []
      if use_multimer:
        model_names += [f"model_{k}_multimer_v3" for k in [1,2,3,4,5]]
      else:
        if use_templates:
          model_names += [f"model_{k}_ptm" for k in [1,2]]
        else:
          model_names += [f"model_{k}_ptm" for k in [1,2,3,4,5]]

    self._model_params, self._model_names = [],[]
    for model_name in model_names:
      params = data.get_model_haiku_params(model_name=model_name, data_dir=data_dir, fuse=True)
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
    def _model(params, model_params, inputs, key):
      inputs["params"] = params
      opt = inputs["opt"]

      aux = {}
      key = Key(key=key).get

      #######################################################################
      # INPUTS
      #######################################################################
      L = inputs["aatype"].shape[0]

      # get sequence
      seq_key = key() if a["shuffle_msa"] else None
      seq = self._get_seq(inputs, aux, seq_key)
            
      # update sequence features      
      pssm = jnp.where(opt["pssm_hard"], seq["hard"], seq["pseudo"])
      if a["use_mlm"]:
        mlm = opt.get("mlm",jax.random.bernoulli(key(),opt["mlm_dropout"],(L,)))
        update_seq(seq["pseudo"], inputs, seq_pssm=pssm, mlm=mlm)
      else:
        update_seq(seq["pseudo"], inputs, seq_pssm=pssm)
      
      # update amino acid sidechain identity
      update_aatype(seq["pseudo"][0].argmax(-1), inputs)      

      inputs["seq"] = aux["seq"]

      # update template features
      inputs["mask_template_interchain"] = opt["template"]["rm_ic"]
      if a["use_templates"]:
        self._update_template(inputs, key())
      
      # set dropout
      inputs["dropout_scale"] = jnp.array(opt["dropout"], dtype=float)

      if "batch" not in inputs:
        inputs["batch"] = None

      # pre callback
      for fn in self._callbacks["model"]["pre"]:
        fn_args = {"inputs":inputs, "opt":opt, "aux":aux,
                   "seq":seq, "key":key(), "params":params}
        sub_args = {k:fn_args.get(k,None) for k in signature(fn).parameters}
        fn(**sub_args)
      
      #######################################################################
      # OUTPUTS
      #######################################################################

      outputs = runner.apply(model_params, key(), inputs)

      # add aux outputs
      aux.update({"atom_positions": outputs["structure_module"]["final_atom_positions"],
                  "atom_mask":      outputs["structure_module"]["final_atom_mask"],                  
                  "residue_index":  inputs["residue_index"],
                  "aatype":         inputs["aatype"],
                  "plddt":          get_plddt(outputs),
                  "pae":            get_pae(outputs), 
                  "ptm":            get_ptm(inputs, outputs),
                  "i_ptm":          get_ptm(inputs, outputs, interface=True), 
                  "cmap":           get_contact_map(outputs, opt["con"]["cutoff"]),
                  "i_cmap":         get_contact_map(outputs, opt["i_con"]["cutoff"]),
                  "prev":           outputs["prev"]})

      #######################################################################
      # LOSS
      #######################################################################
      aux["losses"] = {}

      # add protocol specific losses
      self._get_loss(inputs=inputs, outputs=outputs, aux=aux)

      # sequence entropy loss
      aux["losses"].update(get_seq_ent_loss(inputs))
      
      # experimental masked-language-modeling
      if a["use_mlm"]:
        aux["mlm"] = outputs["masked_msa"]["logits"]
        aux["losses"].update(get_mlm_loss(outputs, mask=mlm, truth=seq["pssm"]))

      # run user defined callbacks
      for c in ["loss","post"]:
        for fn in self._callbacks["model"][c]:
          fn_args = {"inputs":inputs, "outputs":outputs, "opt":opt,
                     "aux":aux, "seq":seq, "key":key(), "params":params}
          sub_args = {k:fn_args.get(k,None) for k in signature(fn).parameters}
          if c == "loss": aux["losses"].update(fn(**sub_args))
          if c == "post": fn(**sub_args)

      # save for debugging
      if a["debug"]: aux["debug"] = {"inputs":inputs,"outputs":outputs}
  
      # weighted loss
      w = opt["weights"]
      loss = sum([v * w[k] if k in w else v for k,v in aux["losses"].items()])
      return loss, aux
    
    return {"grad_fn":jax.jit(jax.value_and_grad(_model, has_aux=True, argnums=0)),
            "fn":jax.jit(_model), "runner":runner}
