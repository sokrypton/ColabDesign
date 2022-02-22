import numpy as np
import jax
import jax.numpy as jnp
from jax.example_libraries.optimizers import sgd

def clear_mem():
  backend = jax.lib.xla_bridge.get_backend()
  for buf in backend.live_buffers(): buf.delete()

from alphafold.common import protein
from alphafold.data import pipeline, templates
from alphafold.model import data, config, model, modules
from alphafold.common import residue_constants

from alphafold.model import all_atom
from alphafold.model import folding

# custom functions
from utils import *
import colabfold as cf

import py3Dmol
import matplotlib
from matplotlib import animation
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

class mk_design_model:
  ######################################
  # model initialization
  ######################################
  def __init__(self, num_seq=1, protocol="fixbb",
               num_models=1, model_mode="sample", model_parallel=False,
               num_recycles=0, recycle_mode="sample",
               use_templates=None):
    
    # decide if templates should be used
    if use_templates is None:
      use_templates = True if protocol=="binder" else False

    self.args = {"num_seq":num_seq, "use_templates":use_templates,
                 "num_models":num_models,
                 "model_mode":model_mode, "model_parallel": model_parallel,
                 "num_recycles":num_recycles, "recycle_mode":recycle_mode}
    
    self.protocol = protocol

    # setup which model params to use
    if use_templates:
      model_name = "model_1_ptm"
      self.args["num_models"] = min(num_models, 2)
    else:
      model_name = "model_3_ptm"
    cfg = config.model_config(model_name)

    # enable checkpointing
    cfg.model.global_config.use_remat = True

    # subbatch_size / chunking
    cfg.model.global_config.subbatch_size = None

    # number of sequences
    if use_templates:
      cfg.data.eval.max_templates = 1
      cfg.data.eval.max_msa_clusters = num_seq + 1
    else:
      cfg.data.eval.max_msa_clusters = num_seq

    cfg.data.common.max_extra_msa = 1
    cfg.data.eval.masked_msa_replace_fraction = 0

    # number of recycles
    cfg.model.num_recycle = num_recycles
    cfg.data.common.num_recycle = num_recycles

    # backprop through recycles
    cfg.model.add_prev = recycle_mode == "add_prev"
    cfg.model.backprop_recycle = recycle_mode == "backprop"
    cfg.model.embeddings_and_evoformer.backprop_dgram = recycle_mode == "backprop"

    self._config = cfg

    # setup model
    self._model_params = [data.get_model_haiku_params(model_name=model_name, data_dir=".")]
    self._runner = model.RunModel(self._config, self._model_params[0], is_training=True)

    # load the other model_params
    if use_templates:
      model_names = ["model_2_ptm"]
    else:
      model_names = ["model_1_ptm","model_2_ptm","model_4_ptm","model_5_ptm"]

    if self.args["model_mode"] == "fixed":
      model_names = model_names[:self.args["num_models"]]

    for model_name in model_names:
      params = data.get_model_haiku_params(model_name, '.')
      self._model_params.append({k: params[k] for k in self._runner.params.keys()})

    # define gradient function
    if self.args["model_parallel"]:
      self._sample_params = jax.jit(lambda y,n:jax.tree_map(lambda x:x[n],y))
      in_axes = (None,0,None,None,None)
      self._model_params = jax.tree_multimap(lambda *x:jnp.stack(x),*self._model_params)
      self._grad_fn, self._fn = [jax.jit(jax.vmap(x,in_axes)) for x in self._get_fn()]

    else:        
      self._grad_fn, self._fn = [jax.jit(x) for x in self._get_fn()]

    # define input function
    if protocol == "fixbb":          self.prep_inputs = self._prep_fixbb
    if protocol == "hallucination":  self.prep_inputs = self._prep_hallucination
    if protocol == "binder":         self.prep_inputs = self._prep_binder

  ######################################
  # setup gradient
  ######################################
  def _get_fn(self):

    # setup function to get gradients
    def mod(params, model_params, inputs, key, opt):

      # initialize the loss function
      losses = {}
      w = opt["weights"]

      # set sequence
      seq_logits = params["seq_logits"] - params["seq_logits"].mean(-1,keepdims=True)

      # shuffle msa
      if self.args["num_seq"] > 1:
        i = jax.random.randint(key,[],0,seq_logits.shape[0])
        seq_logits = seq_logits.at[0].set(seq_logits[i]).at[i].set(seq_logits[0])

      # straight-through/reparameterization
      seq = seq_logits + jnp.where(opt["gumbel"], jax.random.gumbel(key,seq_logits.shape), 0.0)
      seq_soft = jax.nn.softmax(seq / opt["temp"])
      seq_hard = jax.nn.one_hot(seq_soft.argmax(-1),20)
      seq_hard = jax.lax.stop_gradient(seq_hard - seq_soft) + seq_soft

      # create pseudo sequence
      seq_pseudo = opt["soft"] * seq_soft + (1-opt["soft"]) * seq_logits
      seq_pseudo = jnp.where(opt["hard"], seq_hard, seq_pseudo)
      
      # save for aux output
      aux = {"seq":seq_hard,"seq_pseudo":seq_pseudo}

      # entropy loss for msa
      if self.args["num_seq"] > 1:
        seq_prf = seq_hard.mean(0)
        losses["msa_ent"] = -(seq_prf * jnp.log(seq_prf + 1e-8)).sum(-1).mean()
      
      if self.protocol == "binder":
        # concatenate target and binder sequence
        seq_target = jax.nn.one_hot(self._batch["aatype"][:self._target_len],20)
        seq_target = jnp.broadcast_to(seq_target,(self.args["num_seq"],*seq_target.shape))
        seq_pseudo = jnp.concatenate([seq_target, seq_pseudo], 1)
        seq_pseudo = jnp.pad(seq_pseudo,[[0,1],[0,0],[0,0]])
      
      if self.protocol == "hallucination" and self._copies > 1:
        seq_pseudo = jnp.concatenate([seq_pseudo]*self._copies, 1)
      
      # update sequence
      update_seq(seq_pseudo, inputs)
      
      # update amino acid sidechain identity
      N,L = inputs["aatype"].shape[:2]
      aatype = jnp.broadcast_to(jax.nn.one_hot(seq_pseudo[0].argmax(-1),21),(N,L,21))
      update_aatype(aatype, inputs)

      # update template sequence
      if self.protocol == "fixbb" and self.args["use_templates"]:
        # TODO
        inputs["template_aatype"] = inputs["template_aatype"].at[...,:].set(0)
        inputs["template_all_atom_masks"] = inputs["template_all_atom_masks"].at[...,5:].set(0.0)

      # set number of recycles to use
      if self.args["recycle_mode"] == "sample":
        inputs["num_iter_recycling"] = jnp.asarray([opt["recycles"]])
      
      # scale dropout rate
      inputs["scale_rate"] = jnp.where(opt["dropout"],
                                       jnp.full(1,opt["dropout_scale"]),jnp.zeros(1))

      # get outputs
      outputs = self._runner.apply(model_params, key, inputs)
              
      # confidence losses
      pae_prob = jax.nn.softmax(outputs["predicted_aligned_error"]["logits"])
      pae_loss = (pae_prob * jnp.arange(pae_prob.shape[-1])).mean(-1)
      
      plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
      plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)
      
      # define contact
      dgram_bins = jnp.append(0,outputs["distogram"]["bin_edges"])
      if "con_cutoff" in opt:
        con_bins = dgram_bins < opt["con_cutoff"]
      else:
        con_bins = dgram_bins < dgram_bins[-1]
      con_prob = jax.nn.softmax(outputs["distogram"]["logits"])
      con_loss = -jnp.log((con_bins*con_prob).sum(-1))

      # protocol specific losses
      if self.protocol == "binder":

        T = self._target_len
        if self._hotspot is None:
          pae_inter = 0.5 * (pae_loss[...,:T,T:].mean() + pae_loss[...,T:,:T].mean())
          con_inter = con_loss[...,:T,T:].mean()
        else:
          # restrict the inter loss to specific hotspot residues
          H = self._hotspot
          pae_inter = 0.5 * (pae_loss[...,H,T:].mean() + pae_loss[...,T:,H].mean())
          con_inter = con_loss[...,H,T:].mean()

        losses.update({"con_inter":con_inter,"con_intra":con_loss[...,T:,T:].mean(),
                       "pae_inter":pae_inter,"pae_intra":pae_loss[...,T:,T:].mean(),
                       "plddt":plddt_loss[...,T:].mean()})

        aux.update({"pae":get_pae(outputs)})

      if self.protocol == "hallucination":
        losses.update({"con":con_loss.mean(),
                       "plddt":plddt_loss.mean(),
                       "pae":pae_loss.mean()})
        if self._copies > 1:
          aux.update({"pae":get_pae(outputs)})

      if self.protocol == "fixbb":
        fape_loss = get_fape_loss(self._batch, outputs, model_config=self._config)      
        dgram_cce_loss = get_dgram_loss(self._batch, outputs, model_config=self._config)
        losses.update({"plddt":plddt_loss.mean(),"pae":pae_loss.mean(),
                       "dgram_cce":dgram_cce_loss, "fape":fape_loss})
        aux.update({"rmsd":get_rmsd_loss(self._batch, outputs)})

      # loss
      loss = sum([v*w[k] if k in w else v for k,v in losses.items()])

      # save aux outputs
      aux.update({"final_atom_positions":outputs["structure_module"]["final_atom_positions"],
                  "final_atom_mask":outputs["structure_module"]["final_atom_mask"],
                  "plddt":get_plddt(outputs),
                  "losses":losses})

      return loss, (aux)
    
    return jax.value_and_grad(mod, has_aux=True, argnums=0), mod
  
  ######################################
  # input prep functions
  ######################################

  def _prep_features(self, length, template_features=None):
    '''process features'''
    num_seq = self.args["num_seq"]
    sequence = "A" * length
    feature_dict = {
        **pipeline.make_sequence_features(sequence=sequence, description="none",
                                          num_res=length),
        **pipeline.make_msa_features(msas=[length*[sequence]],
                                     deletion_matrices=[num_seq*[[0]*length]])
        }
    if template_features is not None: feature_dict.update(template_features)    
    inputs = self._runner.process_features(feature_dict, random_seed=0)
    if num_seq > 1:
      inputs["msa_row_mask"] = jnp.ones_like(inputs["msa_row_mask"])
      inputs["msa_mask"] = jnp.ones_like(inputs["msa_mask"])
    return inputs

  def _prep_pdb(self, pdb_filename, chain=None):
    '''extract features from pdb'''
    protein_obj = protein.from_pdb_string(pdb_to_string(pdb_filename), chain_id=chain)
    batch = {'aatype': protein_obj.aatype,
             'all_atom_positions': protein_obj.atom_positions,
             'all_atom_mask': protein_obj.atom_mask}

    has_ca = batch["all_atom_mask"][:,0] == 1
    batch = jax.tree_map(lambda x:x[has_ca], batch)
    batch.update(all_atom.atom37_to_frames(**batch))

    template_features = {"template_aatype":jax.nn.one_hot(protein_obj.aatype[has_ca],22)[None],
                         "template_all_atom_masks":protein_obj.atom_mask[has_ca][None],
                         "template_all_atom_positions":protein_obj.atom_positions[has_ca][None],
                         "template_domain_names":np.asarray(["None"])}
    return {"batch":batch,
            "template_features":template_features,
            "residue_index": protein_obj.residue_index[has_ca]}

  def _prep_hotspot(self, hotspot):
    res_idx =  self._inputs["residue_index"][0]
    if len(hotspot) > 0:
      hotspot_set = []
      for idx in hotspot.split(","):
        i,j = idx.split("-") if "-" in idx else (idx,None)
        if j is None: hotspot_set.append(int(i))
        else: hotspot_set += list(range(int(i),int(j)+1))
      hotspot_array = []
      for i in hotspot_set:
        if i in res_idx:
          hotspot_array.append(np.where(res_idx == i)[0][0])
      return jnp.asarray(hotspot_array)
    else:
      return None
  
  # prep functions specific to protocol
  def _prep_binder(self, pdb_filename, chain=None, binder_len=50, hotspot="", **kwargs):
    '''prep inputs for binder design'''
    pdb = self._prep_pdb(pdb_filename, chain=chain)

    target_len = pdb["residue_index"].shape[0]
    self._inputs = self._prep_features(target_len, pdb["template_features"])
    self._inputs["residue_index"][:,:] = pdb["residue_index"]

    self._hotspot = self._prep_hotspot(hotspot)

    total_len = target_len + binder_len
    self._inputs = make_fixed_size(self._inputs, self._runner, total_len)
    self._batch = make_fixed_size(pdb["batch"], self._runner, total_len, batch_axis=False)

    # offset residue index for binder
    self._inputs["residue_index"] = self._inputs["residue_index"].copy()
    self._inputs["residue_index"][:,target_len:] = pdb["residue_index"][-1] + np.arange(binder_len) + 50

    self._inputs["seq_mask"] = np.ones_like(self._inputs["seq_mask"])
    self._inputs["msa_mask"] = np.ones_like(self._inputs["msa_mask"])

    self._target_len = target_len
    self._binder_len = self._len = binder_len
    self._k = -1

    self._default_weights = {"msa_ent":0.01,"plddt":0.1,
                             "pae_intra":0.1,"pae_inter":1.0,
                             "con_intra":0.1,"con_inter":0.5}
    self.restart(**kwargs)

  def _prep_fixbb(self, pdb_filename, chain=None, **kwargs):
    '''prep inputs for fixed backbone design'''
    pdb = self._prep_pdb(pdb_filename, chain=chain)
    length = pdb["residue_index"].shape[0]
    self._inputs = self._prep_features(length, pdb["template_features"])
    # update residue index from pdb
    self._inputs["residue_index"][:,:] = pdb["residue_index"]
    self._batch = pdb["batch"]
    self._len = length
    self._k = -1
    # set weights
    self._default_weights = {"msa_ent":0.01,"dgram_cce":1.0,
                             "fape":0.0,"pae":0.1,"plddt":0.1}
    self.restart(**kwargs)
    
  def _prep_hallucination(self, length=100, copies=1, **kwargs):
    '''prep inputs for hallucination'''
    self._len = length
    self._inputs = self._prep_features(length * copies)
    self._copies = copies
    if copies > 1:
      Ls = [length] * copies
      self._inputs["residue_index"][...,:] = cf.chain_break(self._inputs["residue_index"][0], Ls, 50)
    self._k = -1

    # set weights
    self._default_weights = {"msa_ent":0.01,"pae":1.0,"plddt":1.0,"con":0.5}
    self.restart(**kwargs)

  #################################
  # initialization/restart function
  #################################
  def _setup_optimizer(self, lr=0.1, optimizer=None, **kwargs):
    '''setup which optimizer to use'''
    if optimizer is None: optimizer = sgd
    self._init_fun, self._update_fun, self._get_params = optimizer(lr, **kwargs)
    self._k = 0

  def _init_seq(self, x=None):
    '''initialize sequence'''
    self._key, _key = jax.random.split(self._key)
    shape = (self.args["num_seq"],self._len,20)
    if isinstance(x, np.ndarray) or isinstance(x, jnp.ndarray):
      y = jnp.broadcast_to(x, shape)
    elif isinstance(x, str):
      if x == "gumbel": y = jax.random.gumbel(_key, shape)
      if x == "zeros": y = jnp.zeros(shape)
    else:
      y = 0.01 * jax.random.normal(_key, shape)
    self._params = {"seq_logits":y}
    self._state = self._init_fun(self._params)

  def restart(self, weights=None, seed=None, seq_init=None,
              lr=0.1, optimizer=None, **kwargs):    
    
    # set weights and options
    self.opt = {"weights":self._default_weights.copy()}
    if weights is not None: self.opt["weights"].update(weights)
    self.opt.update({"temp":0.5, "soft":1.0, "hard":True,
                     "dropout":True, "dropout_scale":1.0,
                     "gumbel":False})

    # setup optimizer
    self._setup_optimizer(lr=lr, optimizer=optimizer, **kwargs)    
    self._lr = lr    
    
    # initialize sequence
    if seed is None: seed = np.random.randint(100000)
    self._key = jax.random.PRNGKey(seed)
    self._init_seq(seq_init)

    # initialize trajectory
    self.losses,self._traj = [],{"xyz":[],"seq":[],"plddt":[],"pae":[]}
    self._best_loss, self._best_outs = np.inf, None

  ######################################
  # design function
  ######################################
  def _update_grad(self):
    '''update the gradient'''
    self._key, *_key = jax.random.split(self._key,4)

    # decide number of recycles to use
    if self.args["recycle_mode"] == "sample":
      self.opt["recycles"] = int(jax.random.randint(_key[0],[],0,self.args["num_recycles"]+1))
    else:
      self.opt["recycles"] = self.args["num_recycles"]

    # decide which model params to use
    m = self.args["num_models"]
    ns = jnp.arange(2) if self.args["use_templates"] else jnp.arange(5)
    if self.args["model_mode"] == "fixed" or m == len(ns):
      self._model_num = ns
    elif self.args["model_mode"] == "sample":
      self._model_num = jax.random.choice(_key[1],ns,(m,),replace=False)

    # get current params
    self._params = self._get_params(self._state)
    
    # run in parallel
    if self.args["model_parallel"]:
      p = self._model_params
      if self.args["model_mode"] == "sample":
        p = self._sample_params(p, self._model_num)
      (l,o),g = self._grad_fn(self._params, p, self._inputs, _key[2], self.opt)
      self._grad = jax.tree_map(lambda x: x.mean(0), g)
      self._loss,self._outs = l.mean(),jax.tree_map(lambda x:x[0],o)
      self._losses = jax.tree_map(lambda x: x.mean(0), o["losses"])

    # run in serial
    else:
      _loss, _losses, _outs, _grad = [],[],[],[]
      for n in self._model_num:
        (l,o),g = self._grad_fn(self._params, self._model_params[n],
                                self._inputs, _key[2], self.opt)
        _loss.append(l); _outs.append(o); _grad.append(g)
        _losses.append(o["losses"])

      self._grad = jax.tree_multimap(lambda *v: jnp.mean(jnp.asarray(v), axis=0), *_grad)
      self._loss, self._outs = jnp.mean(jnp.asarray(_loss)), _outs[0]
      self._losses = jax.tree_multimap(lambda *v: jnp.mean(jnp.asarray(v)), *_losses)

  def _save_results(self, save_best=False, verbose=True):
    '''save the results and update trajectory'''

    # save best result
    if save_best and self._loss < self._best_loss:
      self._best_loss = self._loss
      self._best_outs = self._outs
    
    # compile losses
    self._losses.update({"model":self._model_num, "loss":self._loss, **self.opt})
    
    if self.protocol == "fixbb":
      # compute sequence recovery
      seqid = (self._outs["seq"].argmax(-1) == self._batch["aatype"]).mean()
      self._losses.update({"seqid":seqid,"rmsd":self._outs["rmsd"]})

    # save losses
    self.losses.append(self._losses)

    # print losses      
    if verbose:
      I = ["model","recycles"]
      F = ["soft","temp","loss","seqid","msa_ent",
           "pae","pae_intra","pae_inter","plddt",
           "con","con_intra","con_inter","dgram_cce","fape","rmsd"]
      I = " ".join([f"{x}: {self._losses[x]}" for x in I if x in self._losses])
      F = " ".join([f"{x}: {self._losses[x]:.3f}" for x in F if x in self._losses])
      print(f"{self._k}\t{I} {F}")

    if self.args["num_recycles"] == self.opt["recycles"]:
      # save trajectory
      ca_xyz = self._outs["final_atom_positions"][:,1,:]
      traj = {"xyz":ca_xyz,"plddt":self._outs["plddt"],"seq":self._outs["seq_pseudo"]}
      if "pae" in self._outs: traj.update({"pae":self._outs["pae"]})
      for k,v in traj.items(): self._traj[k].append(np.array(v))
      
  def _step(self, weights=None, grad_scale=1.0, **kwargs):
    '''do one step'''
    if weights is not None: self.opt["weights"].update(weights)
    
    # update gradient
    self._update_grad()

    # normalize gradient
    g = self._grad["seq_logits"]
    gn = jnp.linalg.norm(g,axis=(-1,-2),keepdims=True)
    self._grad["seq_logits"] *= grad_scale * jnp.sqrt(self._len)/(gn+1e-7)

    # apply gradient
    self._state = self._update_fun(self._k, self._grad, self._state)
    self._k += 1
    self._save_results(**kwargs)

  ##############################################################################
  # DESIGN FUNCTIONS
  ##############################################################################
  def design(self, iters,
             temp=0.5, e_temp=None,
             soft=False, e_soft=None,
             hard=False, dropout=True, **kwargs):
    self.opt.update({"hard":hard,"dropout":dropout})
    if e_soft is None: e_soft = soft
    if e_temp is None: e_temp = temp
    for i in range(iters):
      self.opt["temp"] = temp * ((e_temp/temp) ** (1/(iters-1))) ** i
      self.opt["soft"] = soft + (e_soft-soft) * i/(iters-1)
      grad_scale = (1-self.opt["soft"]) + self.opt["soft"] * 2 * self.opt["temp"]
      self._step(grad_scale=grad_scale, **kwargs)

  def design_logits(self, iters, **kwargs):
    '''optimize logits'''
    self.design(iters, **kwargs)

  def design_soft(self, iters, **kwargs):
    ''' optimize softmax(logits/temp)'''
    self.design(iters, soft=True, **kwargs)
  
  def design_hard(self, iters, **kwargs):
    ''' optimize argmax(logits)'''
    self.design(iters, soft=True, hard=True, **kwargs)

  def design_2stage(self, soft_iters=100, temp_iters=100, hard_iters=50, **kwargs):
    '''two stage design (soft→hard)'''
    self.design(soft_iters, soft=True, temp=1.0, **kwargs)
    self.design(temp_iters, soft=True, temp=1.0,  e_temp=1e-2, **kwargs)
    self.design(hard_iters, soft=True, temp=1e-2, hard=True, save_best=True, **kwargs)

  def design_3stage(self, soft_iters=300, temp_iters=100, hard_iters=50, **kwargs):
    '''three stage design (logits→soft→hard)'''
    self.design(soft_iters, soft=False, e_soft=True, **kwargs)
    self.design(temp_iters, soft=True, temp=0.5, e_temp=1e-2, **kwargs)
    self.design(hard_iters, soft=True, temp=1e-2, hard=True, save_best=True, **kwargs)
  
  ######################################
  # utils
  ######################################
  def get_seqs(self):
    outs = self._outs if self._best_outs is None else self._best_outs
    x = outs["seq"].argmax(-1)
    return ["".join([order_restype[a] for a in s]) for s in x]
  
  def get_loss(self, x="loss"):
    '''output the loss (for entire trajectory)'''
    return np.array([float(loss[x]) for loss in self.losses])

  def save_pdb(self, filename=None):
    '''save pdb coordinates'''
    outs = self._outs if self._best_outs is None else self._best_outs
    aatype = outs["seq"].argmax(-1)[0]
    if self.protocol == "binder":
      aatype_target = self._batch["aatype"][:self._target_len]
      aatype = np.concatenate([aatype_target,aatype])
    if self.protocol == "hallucination" and self._copies > 1:
      aatype = np.concatenate([aatype] * self._copies)
    p = {"residue_index":self._inputs["residue_index"][0],
         "aatype":aatype,
         "atom_positions":outs["final_atom_positions"],
         "atom_mask":outs["final_atom_mask"]}
    b_factors = outs["plddt"][:,None] * p["atom_mask"]
    p = protein.Protein(**p,b_factors=b_factors)
    pdb_lines = protein.to_pdb(p)
    if filename is None: return pdb_lines
    else:
      with open(filename, 'w') as f: f.write(pdb_lines)
  
  ######################################
  # plotting functions
  ######################################
  def animate(self, s=0, e=None, dpi=100):
    sub_traj = {k:v[s:e] for k,v in self._traj.items()}
    if self.protocol == "fixbb":
      pos_ref = self._batch["all_atom_positions"][:,1,:]
      return make_animation(**sub_traj, pos_ref=pos_ref, dpi=dpi)
    
    if self.protocol == "binder":
      outs = self._outs if self._best_outs is None else self._best_outs
      pos_ref = outs["final_atom_positions"][:,1,:]
      return make_animation(**sub_traj, pos_ref=pos_ref,
                            length=self._target_len, dpi=dpi)    

    if self.protocol == "hallucination":
      outs = self._outs if self._best_outs is None else self._best_outs
      pos_ref = outs["final_atom_positions"][:,1,:]
      if self._copies > 1:
        return make_animation(**sub_traj, pos_ref=pos_ref, length=self._len, dpi=dpi)
      else:
        return make_animation(**sub_traj, pos_ref=pos_ref, dpi=dpi)

  def plot_pdb(self):
    '''use py3Dmol to plot pdb coordinates'''
    view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js')
    view.addModel(self.save_pdb(),'pdb')
    view.setStyle({'cartoon': {}})
    BB = ['C','O','N']
    view.addStyle({'and':[{'resn':["GLY","PRO"],'invert':True},{'atom':BB,'invert':True}]},
                  {'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
    view.addStyle({'and':[{'resn':"GLY"},{'atom':'CA'}]},
                  {'sphere':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
    view.addStyle({'and':[{'resn':"PRO"},{'atom':['C','O'],'invert':True}]},
                  {'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})  
    view.zoomTo()
    view.show()
  
  def plot_traj(self):
    ## TODO
    # currently only fixbb trajories support
    fig, ax1 = plt.subplots(dpi=100)
    ax2 = ax1.twinx()
    if self.protocol == "fixbb":
      ax1.plot(self.get_loss("rmsd"),color="black")
      ax2.plot(self.get_loss("seqid"),color="green",label="seqid")
      # axes labels
      ax1.set_yscale("log")
      ticks = [0.25,0.5,1,2,4,8,16,32,64]
      ax1.set_yticks(ticks);ax1.set_yticklabels(ticks)
      ax1.set_ylabel("rmsd",color="black");ax2.set_ylabel("seqid",color="green")
      ax1.set_xlabel("iterations")
      ax1.set_ylim(0.25,64);ax2.set_ylim(0,0.4)    
      # extras
      ax2.plot(0.4*self.get_loss("soft"),color="yellow",label="soft")
      ax2.plot(0.4*self.get_loss("temp"),color="orange",label="temp")
      ax2.plot(0.4*self.get_loss("hard"),color="red",label="hard")
    #if self.protocol == "hallucination":
    #  ax1.set_ylim(0,1);ax2.set_ylim(0,1)
    #  ax2.plot(self.get_loss("pae"),label="pae")
    #  ax2.plot(self.get_loss("plddt"),label="plddt")

    plt.legend()
    plt.show()

#####################################################################
# UTILS
#####################################################################

def make_animation(xyz, seq, plddt=None, pae=None,
                   pos_ref=None, line_w=2.0,
                   dpi=100, interval=60, color_msa="Taylor",
                   length=None):

  def align(P, Q, P_trim=None):
    if P_trim is None: P_trim = P
    p_trim = P_trim - P_trim.mean(0,keepdims=True)
    p = P - P_trim.mean(0,keepdims=True)
    q = Q - Q.mean(0,keepdims=True)
    return p @ cf.kabsch(p_trim,q)

  # compute reference position
  if pos_ref is None: pos_ref = xyz[-1]
  if length is None: length = len(pos_ref)
  
  # align to reference
  pos_ref_trim = pos_ref[:length]
  # align to reference position
  new_positions = []
  for i in range(len(xyz)):
    new_positions.append(align(xyz[i],pos_ref_trim,xyz[i][:length]))
  pos = np.asarray(new_positions)

  # rotate for best view
  pos_mean = np.concatenate(pos,0)
  m = pos_mean.mean(0)
  rot_mtx = cf.kabsch(pos_mean - m, pos_mean - m, return_v=True)
  pos = (pos - m) @ rot_mtx + m
  pos_ref_full = ((pos_ref - pos_ref_trim.mean(0)) - m) @ rot_mtx + m

  # initialize figure
  if pae is not None and len(pae) == 0: pae = None
  fig = plt.figure()
  gs = GridSpec(4,3, figure=fig)
  if pae is not None:
    ax1, ax2, ax3 = fig.add_subplot(gs[:3,:2]), fig.add_subplot(gs[3:,:]), fig.add_subplot(gs[:3,2:])
  else:
    ax1, ax2 = fig.add_subplot(gs[:3,:]), fig.add_subplot(gs[3:,:])

  fig.subplots_adjust(top=0.95,bottom=0.1,right=0.95,left=0.05,hspace=0,wspace=0)
  fig.set_figwidth(8); fig.set_figheight(6); fig.set_dpi(dpi)
  ax2.set_xlabel("positions"); ax2.set_yticks([])
  if seq[0].shape[0] > 1: ax2.set_ylabel("sequences")
  else: ax2.set_ylabel("amino acids")

  ax1.set_title("N→C") if plddt is None else ax1.set_title("pLDDT")
  if pae is not None:
    ax3.set_title("pAE")
    ax3.set_xticks([])
    ax3.set_yticks([])

  # set bounderies
  x_min,y_min,z_min = np.minimum(np.mean(pos.min(1),0),pos_ref_full.min(0)) - 5
  x_max,y_max,z_max = np.maximum(np.mean(pos.max(1),0),pos_ref_full.max(0)) + 5

  x_pad = ((y_max - y_min) * 2 - (x_max - x_min)) / 2
  y_pad = ((x_max - x_min) / 2 - (y_max - y_min)) / 2
  if x_pad > 0:
    x_min -= x_pad
    x_max += x_pad
  else:
    y_min -= y_pad
    y_max += y_pad

  ax1.set_xlim(x_min, x_max)
  ax1.set_ylim(y_min, y_max)
  ax1.set_xticks([])
  ax1.set_yticks([])

  # get animation frames
  ims = []
  for k in range(len(pos)):
    ims.append([])
    if plddt is None:
      ims[-1].append(cf.plot_pseudo_3D(pos[k], ax=ax1, line_w=line_w, zmin=z_min, zmax=z_max))
    else:
      ims[-1].append(cf.plot_pseudo_3D(pos[k], c=plddt[k], cmin=0.5, cmax=0.9, ax=ax1, line_w=line_w, zmin=z_min, zmax=z_max))
    if seq[k].shape[0] == 1:
      ims[-1].append(ax2.imshow(seq[k][0].T, animated=True, cmap="bwr_r",vmin=-1, vmax=1))
    else:
      cmap = matplotlib.colors.ListedColormap(jalview_color_list[color_msa])
      vmax = len(jalview_color_list[color_msa]) - 1
      ims[-1].append(ax2.imshow(seq[k].argmax(-1), animated=True, cmap=cmap, vmin=0, vmax=vmax, interpolation="none"))
    if pae is not None:
      ims[-1].append(ax3.imshow(pae[k], animated=True, cmap="bwr",vmin=0, vmax=30))

  # make animation!
  ani = animation.ArtistAnimation(fig, ims, blit=True, interval=interval)
  plt.close()
  return ani.to_html5_video()
