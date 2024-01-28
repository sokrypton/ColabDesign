import os
import jax
import jax.numpy as jnp
import numpy as np
from inspect import signature

from colabdesign.af.alphafold.model import data, config, model, all_atom
from colabdesign.af.alphafold.model.msa import make_msa_feats

from colabdesign.shared.model import design_model
from colabdesign.shared.utils import Key

from colabdesign.af.prep   import _af_prep
from colabdesign.af.loss   import _af_loss, get_plddt, get_pae, get_ptm
from colabdesign.af.loss   import get_contact_map
from colabdesign.af.utils  import _af_utils
from colabdesign.af.design import _af_design
from colabdesign.af.inputs import _af_inputs

################################################################
# MK_DESIGN_MODEL - initialize model, and put it all together
################################################################

class mk_af_model(design_model, _af_inputs, _af_loss, _af_prep, _af_design, _af_utils):
  def __init__(self,
               protocol="contigs",
               model_type="alphafold2_ptm",
               use_templates=False,
               debug=False,
               data_dir=".",
               **kwargs):  
    assert protocol in ["fixbb","hallucination","binder","contigs","partial"]
    assert model_type in ["alphafold2","alphafold2_ptm","alphafold2_multimer_v3"]
    if kwargs.pop("use_multimer",False):
      model_type = "alphafold2_multimer_v3"

    self.protocol = protocol
    self._num = kwargs.pop("num_seq",1)
    self._args = {
      # structure options
      "use_templates":use_templates, "num_templates":1, "use_batch_as_template":True,
      "use_initial_guess":False, "use_initial_atom_pos":False,"use_drmsd":False,
      "use_dgram":False, "use_dgram_pred":False, "realign": True, "use_sidechains": False, 
      
      # sequence options
      "optimize_seq":True, "alphabet_size":20, "shuffle_first":True, "copies":1,
      "num_msa":512, "num_extra_msa":1024, "use_mlm": False, "use_cluster_profile": False,

      # model options
      "model_type":model_type,
      "block_diag":not "multimer" in model_type,
      "use_ptm": "ptm" in model_type or "multimer" in model_type:
      "pseudo_multimer":False, 

      # optimizer options
      "optimizer":"sgd", "best_metric":"loss", 
      "traj_iter":1, "traj_max":10000, "clear_prev": True,  "recycle_mode":"last", 

      # technical options
      "use_remat":True, "use_bfloat16":True, "debug":debug
    }

    if self.protocol == "binder": self._args["use_templates"] = True

    self._opt = {"dropout":True, "pssm_hard":False, "learning_rate":0.1, "norm_seq_grad":True,
                 "num_recycles":0, "num_models":1, "sample_models":True, "fape_cutoff":10.0,
                 "temp":1.0, "soft":0.0, "hard":0.0, "alpha":2.0, "partial_loss":True,
                 "con":      {"num":2, "cutoff":14.0, "binary":False, "seqsep":9, "num_pos":float("inf")},
                 "i_con":    {"num":1, "cutoff":21.6875, "binary":False, "num_pos":float("inf")},
                 "template": {"rm":False, "rm_ic":False, "rm_sc":True, "rm_seq":True},                
                 "mlm":      {"replace_fraction":0.15,"uniform_prob":0.1,"profile_prob":0.1,"same_prob":0.1}}
    
    # set arguments/options
    if "initial_guess" in kwargs: kwargs["use_initial_guess"] = kwargs.pop("initial_guess")
    model_names = kwargs.pop("model_names",None)
    keys = list(kwargs.keys())
    for k in keys:
      if k in self._args: self._args[k] = kwargs.pop(k)
      if k in self._opt:  self._opt[k]  = kwargs.pop(k)

    # set weights
    weights_keys = ["plddt","exp_res","helix","con","i_con",
                    "rmsd","dgram_cce", "fape", "sc_fape","sc_chi","sc_chi_norm"]
    if self._args["use_ptm"]: weights_keys += ["pae","i_pae"]
    self._opt["weights"] = {k:0.0 for k in weights_keys}

    self._params = {}
    self._inputs = {}
    self._tmp = {"traj":{"seq":[],"xyz":[],"plddt":[]},"log":[],"best":{}}
    if self._args["use_ptm"]: self._tmp["traj"]["pae"] = []

    self.opt = self._opt # backward compatibility

    assert self._args["recycle_mode"] in ["average","first","last","sample"]

    # collect callbacks
    self._callbacks = {"model": {"pre": kwargs.pop("pre_callback",None),
                                 "post":kwargs.pop("post_callback",None),
                                 "loss":kwargs.pop("loss_callback",None)},
                       "design":{"pre": kwargs.pop("pre_design_callback",None),
                                 "post":kwargs.pop("post_design_callback",None)}}
    
    for m,n in self._callbacks.items():
      for k,v in n.items():
        if v is None: v = []
        if not isinstance(v,list): v = [v]
        self._callbacks[m][k] = v

    old_params = kwargs.pop("old_params",{})
    assert len(kwargs) == 0, f"ERROR: the following inputs were not set: {kwargs}"

    #############################
    # configure AlphaFold
    #############################
    if "multimer" in model_type:
      self._cfg = config.model_config("model_1_multimer")   
      self._cfg.model.embeddings_and_evoformer.pseudo_multimer = self._args["pseudo_multimer"]
      self.opt["pssm_hard"] = True
    else:
      if "ptm" in model_type:
        self._cfg = config.model_config("model_1_ptm" if self._args["use_templates"] else "model_3_ptm")
      else:
        self._cfg = config.model_config("model_1" if self._args["use_templates"] else "model_3")

    self._cfg.model.num_recycle = 0
    self._cfg.model.global_config.use_remat = self._args["use_remat"]
    self._cfg.model.global_config.use_dgram_pred = self._args["use_dgram_pred"]
    self._cfg.model.global_config.bfloat16  = self._args["use_bfloat16"]
    self._cfg.model.global_config.bfloat16_output  = self._args["use_bfloat16"]

    # load model_params
    if model_names is None:
      model_names = []
      if "multimer" in model_type:
        model_names += [f"model_{k}_multimer_v3" for k in [1,2,3,4,5]]
      elif "ptm" in model_type:
        if self._args["use_templates"]:
          model_names += [f"model_{k}_ptm" for k in [1,2]]
        else:
          model_names += [f"model_{k}_ptm" for k in [1,2,3,4,5]]
      else:
        if self._args["use_templates"]:
          model_names += [f"model_{k}" for k in [1,2]]
        else:
          model_names += [f"model_{k}" for k in [1,2,3,4,5]]


    self._model_params, self._model_names = [],[]
    for model_name in model_names:
      if model_name in old_params:
        params = old_params[model_name]
      else:
        params = data.get_model_haiku_params(model_name=model_name, data_dir=data_dir,
          fuse=True, rm_templates=(not self._args["use_templates"] and not "multimer" in model_type))
      if params is not None:
        self._model_params.append(params)
        self._model_names.append(model_name)
      else:
        print(f"WARNING: '{model_name}' not found")

    #####################################
    # set protocol specific functions
    #####################################
    idx = ["fixbb","hallucination","binder","contigs","partial"].index(self.protocol)
    self.prep_inputs = [self._prep_fixbb, self._prep_hallucination, self._prep_binder, self._prep_contigs, self._prep_partial][idx]

  def _get_model(self, cfg, callback=None):

    a = self._args
    runner = model.RunModel(cfg, use_multimer="multimer" in model_type)

    # setup function to get gradients
    def _model(params, model_params, inputs, key):
            
      opt = inputs["opt"]
      aux = {}
      key = Key(key=key).get

      #######################################################################
      # INPUTS
      #######################################################################
      # get sequence
      use_cluster_profile = (not self._args["optimize_seq"] and self._args["use_cluster_profile"])
      self._update_seq(params, inputs, aux, key())
      inputs = make_msa_feats(inputs, key(),
        num_msa=self._args["num_msa"], num_extra_msa=self._args["num_extra_msa"],
        use_mlm=self._args["use_mlm"], mlm_opt=opt["mlm"],
        use_cluster_profile=use_cluster_profile)
      seq = aux["seq"] = inputs["msa"]

      # update template features
      inputs["mask_template_interchain"] = opt["template"]["rm_ic"]
      if a["use_templates"]:
        self._update_template(inputs)
      
      # set dropout
      inputs["use_dropout"] = opt["dropout"]

      if "batch" not in inputs:
        inputs["batch"] = None

      # optimize model params
      if "model_params" in params:
        model_params.update(params["model_params"])

      # pre callback
      for fn in self._callbacks["model"]["pre"]:
        fn_args = {"inputs":inputs, "opt":opt, "aux":aux, "seq":seq, 
                   "key":key(), "params":params, "model_params":model_params}
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
                  "cmap":           get_contact_map(outputs, opt["con"]["cutoff"]),
                  "i_cmap":         get_contact_map(outputs, opt["i_con"]["cutoff"]),
                  "prev":           outputs["prev"]})

      if a["use_ptm"]:
        aux.update({
          "pae":          get_pae(outputs),
          "ptm":          get_ptm(inputs, outputs),
          "i_ptm":        get_ptm(inputs, outputs, interface=True)
        })

      #######################################################################
      # LOSS
      #######################################################################
      aux["losses"] = {}

      # add protocol specific losses
      self._get_loss(inputs=inputs, outputs=outputs, aux=aux)

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
