import random, copy, os
import numpy as np
import jax
import jax.numpy as jnp

try:
  from jax.example_libraries.optimizers import sgd, adam
except:
  from jax.experimental.optimizers import sgd, adam

def clear_mem():
  backend = jax.lib.xla_bridge.get_backend()
  for buf in backend.live_buffers(): buf.delete()

from alphafold.common import protein
from alphafold.data import pipeline, templates
from alphafold.model import data, config, model, modules
from alphafold.common import residue_constants

from alphafold.model import all_atom
from alphafold.model import folding
from alphafold.data import prep_inputs

# custom functions
from utils import *
import colabfold as cf

import py3Dmol
import matplotlib
from matplotlib import animation
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
 
#################################################
# AF_INIT - input prep functions
#################################################
class _af_init:
  def repeat_idx(self, idx, copies=1, offset=50):
    idx_offset = np.repeat(np.cumsum([0]+[idx[-1]+offset]*(copies-1)),len(idx))
    return np.tile(idx,copies) + idx_offset

  def _prep_features(self, length, template_features=None):
    '''process features'''
    num_seq = self.args["num_seq"]
    sequence = "A" * length
    feature_dict = {
        **pipeline.make_sequence_features(sequence=sequence, description="none", num_res=length),
        **pipeline.make_msa_features(msas=[length*[sequence]], deletion_matrices=[num_seq*[[0]*length]])
    }
    if template_features is not None: feature_dict.update(template_features)    
    inputs = self._runner.process_features(feature_dict, random_seed=0)
    if num_seq > 1:
      inputs["msa_row_mask"] = jnp.ones_like(inputs["msa_row_mask"])
      inputs["msa_mask"] = jnp.ones_like(inputs["msa_mask"])
    return inputs

  def _prep_pdb(self, pdb_filename, chain=None):
    '''extract features from pdb'''
    if chain is None: chains = [None]
    else: chains = chain.split(",")
    o,last = [],0
    for chain in chains:
      protein_obj = protein.from_pdb_string(pdb_to_string(pdb_filename), chain_id=chain)
      batch = {'aatype': protein_obj.aatype,
               'all_atom_positions': protein_obj.atom_positions,
               'all_atom_mask': protein_obj.atom_mask}

      has_ca = batch["all_atom_mask"][:,0] == 1
      batch = jax.tree_map(lambda x:x[has_ca], batch)
      batch.update(all_atom.atom37_to_frames(**batch))

      seq = "".join([order_restype[a] for a in batch["aatype"]])
      template_aatype = residue_constants.sequence_to_onehot(seq, residue_constants.HHBLITS_AA_TO_ID)
      template_features = {"template_aatype":template_aatype,
                           "template_all_atom_masks":batch["all_atom_mask"],
                           "template_all_atom_positions":batch["all_atom_positions"]}
      
      residue_index = protein_obj.residue_index[has_ca] + last
      last = residue_index[-1] + 50
      o.append({"batch":batch,
                "template_features":template_features,
                "residue_index": residue_index})
    o = jax.tree_multimap(lambda *x:np.concatenate(x,0),*o)
    o["template_features"] = {"template_domain_names":np.asarray(["None"]),
                              **jax.tree_map(lambda x:x[None],o["template_features"])}
    return o                          

  def _prep_pos(self, pos, res_idx=None):
    if res_idx is None:
      res_idx = self._inputs["residue_index"][0]
    if len(pos) > 0:
      pos_set = []
      for idx in pos.split(","):
        i,j = idx.split("-") if "-" in idx else (idx,None)
        if j is None: pos_set.append(int(i))
        else: pos_set += list(range(int(i),int(j)+1))
      pos_array = []
      for i in pos_set:
        if i in res_idx:
          pos_array.append(np.where(res_idx == i)[0][0])
      return jnp.asarray(pos_array)
    else:
      return None
  
  # prep functions specific to protocol
  def _prep_binder(self, pdb_filename, chain=None, binder_len=50, hotspot="", **kwargs):
    '''prep inputs for binder design'''

    # get pdb info
    pdb = self._prep_pdb(pdb_filename, chain=chain)
    target_len = pdb["residue_index"].shape[0]
    self._inputs = self._prep_features(target_len, pdb["template_features"])
    self._inputs["residue_index"][...,:] = pdb["residue_index"]

    # gather hotspot info
    self._hotspot = self._prep_pos(hotspot)

    # pad inputs
    total_len = target_len + binder_len
    self._inputs = make_fixed_size(self._inputs, self._runner, total_len)
    self._batch = make_fixed_size(pdb["batch"], self._runner, total_len, batch_axis=False)

    # offset residue index for binder
    self._inputs["residue_index"] = self._inputs["residue_index"].copy()
    self._inputs["residue_index"][:,target_len:] = pdb["residue_index"][-1] + np.arange(binder_len) + 50
    for k in ["seq_mask","msa_mask"]: self._inputs[k] = np.ones_like(self._inputs[k])

    self._target_len = target_len
    self._binder_len = self._len = binder_len

    self._default_weights.update({"con":0.5, "i_pae":0.01, "i_con":0.5})
    self.restart(set_defaults=True, **kwargs)

  def _prep_fixbb(self, pdb_filename, chain=None, copies=1, homooligomer=False, **kwargs):
    '''prep inputs for fixed backbone design'''
    pdb = self._prep_pdb(pdb_filename, chain=chain)
    self._batch = pdb["batch"]
    self._wt_aatype = self._batch["aatype"]
    self._len = pdb["residue_index"].shape[0]
    self._inputs = self._prep_features(self._len, pdb["template_features"])
    self._copies = copies
    
    # set weights
    self._default_weights.update({"dgram_cce":1.0, "fape":0.0, "rmsd":0.0, "con":0.0})

    # update residue index from pdb
    if copies > 1:
      if homooligomer:
        self._len = self._len // copies
        self._inputs["residue_index"] = pdb["residue_index"][None]
      else:
        self._inputs = make_fixed_size(self._inputs, self._runner, self._len * copies)
        self._batch = make_fixed_size(self._batch, self._runner, self._len * copies, batch_axis=False)
        self._inputs["residue_index"] = self.repeat_idx(pdb["residue_index"], copies)[None]
        for k in ["seq_mask","msa_mask"]: self._inputs[k] = np.ones_like(self._inputs[k])
        
      self._default_weights.update({"i_pae":0.01, "i_con":0.0})
    else:
      self._inputs["residue_index"] = pdb["residue_index"][None]

    self.restart(set_defaults=True, **kwargs)
    
  def _prep_hallucination(self, length=100, copies=1, **kwargs):
    '''prep inputs for hallucination'''
    self._len = length
    self._inputs = self._prep_features(length * copies)
    self._copies = copies
    # set weights
    self._default_weights.update({"con":1.0})
    if copies > 1:
      self._inputs["residue_index"] = self.repeat_idx(np.arange(length), copies)[None]
      self._default_weights.update({"i_pae":0.01, "i_con":0.1})

    self.restart(set_defaults=True, **kwargs)

  def _prep_partial(self, pdb_filename, chain=None, pos=None, length=None,
                    fix_seq=True, mainchain=True, sidechain=False, **kwargs):
    '''prep input for partial hallucination'''
    self.args.update({"mainchain":mainchain, "sidechain":sidechain, "fix_seq":fix_seq})
    self._copies = 1
    
    # get [pos]itions of interests
    pdb = self._prep_pdb(pdb_filename, chain=chain)
    if pos is None:
      self._default_opt["pos"] = jnp.arange(pdb["residue_index"].shape[0])
    else:
      self._default_opt["pos"] = self._prep_pos(pos, pdb["residue_index"])

    pos = np.asarray(self._default_opt["pos"])
    self._batch = jax.tree_map(lambda x:x[pos], pdb["batch"])
    self._wt_aatype = self._batch["aatype"]

    if sidechain:
      self._batch.update(prep_inputs.make_atom14_positions(self._batch))

    self._len = pdb["residue_index"].shape[0] if length is None else length
    if self.args["use_templates"]:
      temp_features = {'template_aatype': np.zeros([1,self._len,22]),
                       'template_all_atom_masks': np.zeros([1,self._len,37]),
                       'template_all_atom_positions': np.zeros([1,self._len,37,3]),
                       'template_domain_names': np.asarray(["None"])}
      self._inputs = self._prep_features(self._len, temp_features)
      if fix_seq:
        self._default_opt["template_aatype"] = self._batch["aatype"]

    else:
      self._inputs = self._prep_features(self._len)
    
    weights = {"dgram_cce":1.0,"con":1.0}
    if mainchain: weights.update({"fape":0.0, "rmsd":0.0})
    if sidechain: weights.update({"sc_fape":0.0, "sc_rmsd":0.0})
    self._default_weights.update(weights)
    
    self.restart(set_defaults=True, **kwargs)

####################################################
# AF_LOSS - setup loss function
####################################################
class _af_loss:
  def _get_fn(self):
    # setup function to get gradients
    def mod(params, model_params, inputs, key, opt):

      # initialize the loss function
      losses = {}

      # set sequence
      seq = params["seq"]

      # shuffle msa (randomly pick which sequence is query)
      if self.args["num_seq"] > 1:
        i = jax.random.randint(key,[],0,seq.shape[0])
        seq = seq.at[0].set(seq[i]).at[i].set(seq[0])

      # straight-through/reparameterization
      seq_logits = 2.0 * seq + opt["bias"] + jnp.where(opt["gumbel"], jax.random.gumbel(key,seq.shape), 0.0)
      seq_soft = jax.nn.softmax(seq_logits / opt["temp"])
      seq_hard = jax.nn.one_hot(seq_soft.argmax(-1), 20)
      seq_hard = jax.lax.stop_gradient(seq_hard - seq_soft) + seq_soft

      # create pseudo sequence
      seq_pseudo = opt["soft"] * seq_soft + (1-opt["soft"]) * seq
      seq_pseudo = opt["hard"] * seq_hard + (1-opt["hard"]) * seq_pseudo

      # for partial hallucination, force sequence if sidechain constraints defined
      if self.protocol == "partial" and (self.args["fix_seq"] or self.args["sidechain"]):
        seq_ref = jax.nn.one_hot(self._wt_aatype,20)
        seq_pseudo = seq_pseudo.at[:,opt["pos"],:].set(seq_ref)
        seq_hard = seq_hard.at[:,opt["pos"],:].set(seq_ref)

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
        if self.args["use_templates"]:
          seq_pseudo = jnp.pad(seq_pseudo,[[0,1],[0,0],[0,0]])
      
      if self.protocol in ["fixbb","hallucination"] and self._copies > 1:
        seq_pseudo = jnp.concatenate([seq_pseudo]*self._copies, 1)
            
      # update sequence
      update_seq(seq_pseudo, inputs)
      
      # update amino acid sidechain identity
      N,L = inputs["aatype"].shape[:2]
      aatype = jnp.broadcast_to(jax.nn.one_hot(seq_pseudo[0].argmax(-1),21),(N,L,21))
      update_aatype(aatype, inputs)

      # update template features
      if self.args["use_templates"]:

        if self.protocol == "partial":
          pb, pb_mask = model.modules.pseudo_beta_fn(self._batch["aatype"],
                                                     self._batch["all_atom_positions"],
                                                     self._batch["all_atom_mask"])       
          
          template_feats = {"template_aatype":opt["template_aatype"],
                            "template_all_atom_positions":self._batch["all_atom_positions"],
                            "template_all_atom_masks":self._batch["all_atom_mask"],
                            "template_pseudo_beta":pb,
                            "template_pseudo_beta_mask":pb_mask}

          for k,v in template_feats.items():
            inputs[k] = inputs[k].at[:,:,opt["pos"]].set(v)

        if self.protocol == "fixbb":
          inputs["template_aatype"] = inputs["template_aatype"].at[:].set(opt["template_aatype"])
        
        # disable sidechains
        inputs["template_all_atom_masks"] = inputs["template_all_atom_masks"].at[...,5:].set(0.0)

        # dropout template input features
        key, key_ = jax.random.split(key)
        pos_mask = jax.random.bernoulli(key_,1-opt["template_dropout"],(L,))
        inputs["template_all_atom_masks"] = inputs["template_all_atom_masks"] * pos_mask[:,None]
        inputs["template_pseudo_beta_mask"] = inputs["template_pseudo_beta_mask"] * pos_mask

      # set number of recycles to use
      if self.args["recycle_mode"] not in ["backprop","add_prev"]:
        inputs["num_iter_recycling"] = jnp.asarray([opt["recycles"]])
      
      # scale dropout rate
      inputs["scale_rate"] = jnp.where(opt["dropout"],jnp.full(1,opt["dropout_scale"]),jnp.zeros(1))
      
      # get outputs
      outputs = self._runner.apply(model_params, key, inputs)

      if self.args["recycle_mode"] == "average":
        aux["init"] = {'init_pos': outputs['structure_module']['final_atom_positions'][None],
                       'init_msa_first_row': outputs['representations']['msa_first_row'][None],
                       'init_pair': outputs['representations']['pair'][None]}
                    
      # pairwise loss
      o = self._get_pairwise_loss(outputs, opt)
      losses.update(o["losses"]); aux.update(o["aux"])

      # confidence losses      
      plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
      plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)

      if self.protocol == "binder":
        losses["plddt"] = plddt_loss[...,self._target_len:].mean()
      else:
        losses["plddt"] = plddt_loss.mean()

      if self.protocol == "fixbb":
        o = self._get_fixbb_loss(outputs, opt)
        losses.update(o["losses"]); aux.update(o["aux"])

      if self.protocol == "partial":
        o = self._get_partial_loss(outputs, opt)
        losses.update(o["losses"]); aux.update(o["aux"])

      # weighted loss
      loss = sum([v*opt["weights"][k] for k,v in losses.items()])

      # save aux outputs
      aux.update({"final_atom_positions":outputs["structure_module"]["final_atom_positions"],
                  "final_atom_mask":outputs["structure_module"]["final_atom_mask"],
                  "plddt":get_plddt(outputs), "losses":losses})

      return loss, (aux)
    
    return jax.value_and_grad(mod, has_aux=True, argnums=0), mod

  def _get_fixbb_loss(self, outputs, outs):
    losses, aux = {},{}
    fape_loss = get_fape_loss(self._batch, outputs, model_config=self._config)      
    dgram_cce_loss = get_dgram_loss(self._batch, outputs, model_config=self._config, copies=self._copies)
    rmsd_loss = get_rmsd_loss_w(self._batch, outputs, copies=self._copies)
    losses.update({"dgram_cce":dgram_cce_loss,
                   "fape":fape_loss,
                   "rmsd":rmsd_loss})
    return {"losses":losses, "aux":aux}    

  def _get_partial_loss(self, outputs, opt):
    '''compute loss for partial hallucination protocol'''
    losses, aux = {},{}
    pos = opt["pos"]
    _config = self._config.model.heads.structure_module

    # dgram
    dgram = outputs["distogram"]["logits"][:,pos][pos,:]
    losses["dgram_cce"] = get_dgram_loss(self._batch, outputs, model_config=self._config, logits=dgram)

    # rmsd
    ca_true = self._batch["all_atom_positions"][:,1]
    ca_pred = outputs["structure_module"]["final_atom_positions"][pos,1]
    losses["rmsd"] = jnp_rmsd(ca_true, ca_pred)
    
    # fape
    fape_loss = {"loss":0.0}
    struct = outputs["structure_module"]
    struct["traj"] = struct["traj"][...,pos,:]
    folding.backbone_loss(fape_loss, self._batch, struct, _config)
    losses["fape"] = fape_loss["loss"]

    # sidechain specific losses
    if self.args["sidechain"]:

      pred_pos = struct["final_atom14_positions"][pos]

      # sc_fape
      struct.update(folding.compute_renamed_ground_truth(self._batch, pred_pos))
      for k in ["frames","atom_pos"]:
        struct["sidechains"][k] = jax.tree_map(lambda x:x[:,pos], struct["sidechains"][k])
      losses["sc_fape"] = folding.sidechain_loss(self._batch, struct, _config)["loss"]

      # sc_rmsd
      true_aa = self._batch["aatype"]
      true_pos = all_atom.atom37_to_atom14(self._batch["all_atom_positions"],self._batch)
      losses["sc_rmsd"] = get_sc_rmsd(true_pos, pred_pos, true_aa)

    return {"losses":losses, "aux":aux}

  def _get_pairwise_loss(self, outputs, opt):
    '''compute pairwise losses (contact and pae)'''

    losses, aux = {},{}
    # pae loss
    pae_prob = jax.nn.softmax(outputs["predicted_aligned_error"]["logits"])
    pae = (pae_prob * jnp.arange(pae_prob.shape[-1])).mean(-1)
    
    # define distogram
    dgram = outputs["distogram"]["logits"]
    dgram_bins = jnp.append(0,outputs["distogram"]["bin_edges"])
    aux.update({"dgram":dgram})

    # contact
    c,ic = opt["con"],opt["i_con"]
    def get_con(x, cutoff, binary=True):
      con_bins = dgram_bins < cutoff
      con_prob = jax.nn.softmax(x - 1e7 * (1-con_bins))
      con_loss_cce = -(con_prob * jax.nn.log_softmax(x)).sum(-1)
      con_loss_bce = -jnp.log((con_bins * jax.nn.softmax(x) + 1e-8).sum(-1))
      return jnp.where(binary, con_loss_bce, con_loss_cce)

    def set_diag(x, k, val=0.0):
      L = x.shape[-1]
      mask = jnp.abs(jnp.arange(L)[:,None] - jnp.arange(L)[None,:]) < k
      return jnp.where(mask, val, x)      
    
    def min_k(x, k=1):
      x = jnp.sort(x)
      mask = jnp.arange(x.shape[-1]) < k
      return jnp.where(mask, x, 0.0).sum(-1)/mask.sum(-1)
    
    # if more than 1 chain, split pae/con into inter/intra
    if self.protocol == "binder" or self._copies > 1:
      aux["pae"] = get_pae(outputs)

      L = self._target_len if self.protocol == "binder" else self._len
      H = self._hotspot if hasattr(self,"_hotspot") else None

      for k,v in zip(["pae","con"],[pae,dgram]):
        
        # split pairwise features into intra (x) and inter (ix)
        aa, bb = v[:L,:L], v[L:,L:]
        ab = v[:L,L:] if H is None else v[H,L:]
        ba = v[L:,:L] if H is None else v[L:,H]
        abba = (ab + ba.swapaxes(0,1)) / 2

        if self.protocol == "binder":
          x,ix = bb,abba.swapaxes(0,1)
        else:
          x,ix = aa,abba
        
        if k == "con":
          losses["helix"] = jnp.diagonal(get_con(x,6.0),3).mean()
          x = get_con(x,c["cutoff"],c["binary"])
          x = min_k(set_diag(x,c["seqsep"],1e8),c["num"])
          ix = get_con(ix,ic["cutoff"],ic["binary"])
          ix = min_k(ix,ic["num"])

        losses.update({f"i_{k}":ix.mean(), k:x.mean()})
    
    else:
      x = get_con(dgram,c["cutoff"],c["binary"])
      x = min_k(set_diag(x,c["seqsep"],1e8),c["num"])
      losses.update({"con":x.mean(),"pae":pae.mean(),
                     "helix":jnp.diagonal(get_con(dgram,8.0),3).mean()})
    
    return {"losses":losses, "aux":aux}

####################################################
# AF_DESIGN - design functions
####################################################
class _af_design:
  def _setup_optimizer(self, optimizer="sgd", lr_scale=1.0, **kwargs):
    '''setup which optimizer to use'''
    if optimizer == "adam":
      optimizer = adam
      lr = 0.02 * lr_scale
    elif optimizer == "sgd":
      optimizer = sgd
      lr = 0.1 * lr_scale
    else:
      lr = lr_scale

    self._init_fun, self._update_fun, self._get_params = optimizer(lr)
    self._k = 0

  def _init_seq(self, mode=None, seq=None, rm_aa=None, add_seq=False, **kwargs):
    '''initialize sequence'''

    # backward compatibility
    if "seq_init" in kwargs:
      x = kwargs["seq_init"]
      if isinstance(x, np.ndarray) or isinstance(x, jnp.ndarray): seq = x
      elif isinstance(x, str):
        if len(x) == self._len: seq = x
        else: mode = x

    if mode is None: mode = ""
    
    # initialize sequence
    self._key, _key = jax.random.split(self._key)
    shape = (self.args["num_seq"], self._len, 20)
    x = 0.01 * jax.random.normal(_key, shape)

    # initialize bias
    b = jnp.zeros(shape[1:])

    if "gumbel" in mode:
      x = jax.random.gumbel(_key, shape) / 2.0

    if "soft" in mode:
      x = jax.nn.softmax(x * 2.0)

    if "wildtype" in mode or "wt" in mode:
      wt = jax.nn.one_hot(self._wt_aatype, 20)
      if self.protocol == "fixbb":
        wt = wt[...,:self._len]
      if self.protocol == "partial":
        x = x.at[...,self.opt["pos"],:].set(wt)
        if add_seq: b = b.at[self.opt["pos"],:].set(wt * 1e8)
      else:
        x = x.at[...,:,:].set(wt)
        if add_seq: b = b.at[:,:].set(wt * 1e8)

    if seq is not None:
      if isinstance(seq, np.ndarray) or isinstance(seq, jnp.ndarray):
        x_ = seq
      if isinstance(seq, str):
        x_ = jnp.array([residue_constants.restype_order.get(aa,-1) for aa in seq])
        x_ = jax.nn.one_hot(x_,20)
      x = x + x_
      if add_seq: b = b + x_ * 1e8

    if rm_aa is not None:
      for aa in rm_aa.split(","):
        b = b - 1e8 * np.eye(20)[residue_constants.restype_order[aa]]

    self.opt["bias"] = b
    self._params = {"seq":x}
    self._state = self._init_fun(self._params)

  def restart(self, seed=None, weights=None, opt=None, set_defaults=False, **kwargs):
    
    # set weights and options
    if set_defaults:
      if weights is not None: self._default_weights.update(weights)
      if opt is not None: self._default_opt.update(opt)

    self.opt = copy.deepcopy(self._default_opt)
    self.opt["weights"] = self._default_weights.copy()

    if not set_defaults:
      if weights is not None: self.opt["weights"].update(weights)
      if opt is not None: self.opt.update(opt)

    # setup optimizer
    self._setup_optimizer(**kwargs)
    
    # initialize sequence
    self._seed = random.randint(0,2147483647) if seed is None else seed
    self._key = jax.random.PRNGKey(self._seed)
    self._init_seq(**kwargs)

    # initialize trajectory
    self.losses,self._traj = [],{"xyz":[],"seq":[],"plddt":[],"pae":[]}
    self._best_loss, self._best_outs = np.inf, None
    
  def _recycle(self, model_params, key, params=None, opt=None):
    if params is None: params = self._params
    if opt is None: opt = self.opt
    if self.args["recycle_mode"] == "average":
      #----------------------------------------
      # average gradients across all recycles
      #----------------------------------------
      L = self._inputs["residue_index"].shape[-1]
      self._inputs.update({'init_pos': np.zeros([1, L, 37, 3]),
                           'init_msa_first_row': np.zeros([1, L, 256]),
                           'init_pair': np.zeros([1, L, L, 128])})
      
      grad = []
      for _ in range(opt["recycles"]+1):
        key, _key = jax.random.split(key)
        (loss,outs),_grad = self._grad_fn(params, model_params, self._inputs, _key, opt)
        grad.append(_grad)
        self._inputs.update(outs["init"])
      grad = jax.tree_multimap(lambda *x: jnp.asarray(x).mean(0), *grad)
      return (loss, outs), grad
    else:
      #----------------------------------------
      # use gradients from last recycle
      #----------------------------------------
      if self.args["use_init_pos"]:
        self._inputs["init_pos"] = self._batch["all_atom_positions"][None]        
      if self.args["recycle_mode"] == "sample":
        # decide number of recycles to use
        key, _key = jax.random.split(key)
        opt["recycles"] = int(jax.random.randint(_key,[],0,self._default_opt["recycles"]+1))
      
      return self._grad_fn(params, model_params, self._inputs, key, opt)
      
  def _update_grad(self):
    '''update the gradient'''

    # get current params
    self._params = self._get_params(self._state)

    # update key
    self._key, key = jax.random.split(self._key)

    # decide which model params to use
    m = self.opt["models"]
    ns = jnp.arange(2) if self.args["use_templates"] else jnp.arange(5)
    if self.args["model_mode"] == "fixed" or m == len(ns):
      self._model_num = ns[:m]
    elif self.args["model_mode"] == "sample":
      key,_key = jax.random.split(key)
      self._model_num = jax.random.choice(_key,ns,(m,),replace=False)
    
    #------------------------
    # run in parallel
    #------------------------
    if self.args["model_parallel"]:
      p = self._model_params
      if self.args["model_mode"] == "sample":
        p = self._sample_params(p, self._model_num)
      (l,o),g = self._recycle(p, key)
      self._grad = jax.tree_map(lambda x: x.mean(0), g)
      self._loss, self._outs = l.mean(),jax.tree_map(lambda x:x[0],o)
      self._losses = jax.tree_map(lambda x: x.mean(0), o["losses"])

    #------------------------
    # run in serial
    #------------------------
    else:
      self._model_num = np.array(self._model_num).tolist()    
      _loss, _losses, _outs, _grad = [],[],[],[]
      for n in self._model_num:
        p = self._model_params[n]
        (l,o),g = self._recycle(p, key)
        _loss.append(l); _outs.append(o); _grad.append(g)
        _losses.append(o["losses"])
      self._grad = jax.tree_multimap(lambda *v: jnp.asarray(v).mean(0), *_grad)      
      self._loss, self._outs = jnp.mean(jnp.asarray(_loss)), _outs[0]
      self._losses = jax.tree_multimap(lambda *v: jnp.asarray(v).mean(), *_losses)

  #-------------------------------------
  # STEP FUNCTION
  #-------------------------------------
  def _step(self, weights=None, lr_scale=1.0, **kwargs):
    '''do one step'''

    if weights is not None: self.opt["weights"].update(weights)
    
    # update gradient
    self._update_grad()

    # normalize gradient
    g = self._grad["seq"]
    gn = jnp.linalg.norm(g,axis=(-1,-2),keepdims=True)
    self._grad["seq"] *= lr_scale * jnp.sqrt(self._len)/(gn+1e-7)

    # apply gradient
    self._state = self._update_fun(self._k, self._grad, self._state)
    self._k += 1
    self._save_results(**kwargs)


  #-----------------------------------------------------------------------------
  # DESIGN FUNCTIONS
  #-----------------------------------------------------------------------------
  def design(self, iters,
              temp=1.0, e_temp=None,
              soft=False, e_soft=None,
              hard=False, dropout=True, gumbel=False, **kwargs):

    self.opt.update({"hard":1.0 if hard else 0.0,"dropout":dropout,"gumbel":gumbel})
    if e_soft is None: e_soft = soft
    if e_temp is None: e_temp = temp
    for i in range(iters):
      self.opt["temp"] = e_temp + (temp - e_temp) * (1-i/(iters-1)) ** 2
      self.opt["soft"] = soft + (e_soft - soft) * i/(iters-1)
      # decay learning rate based on temperature
      lr_scale = (1 - self.opt["soft"]) + (self.opt["soft"] * self.opt["temp"])
      self._step(lr_scale=lr_scale, **kwargs)

  def design_logits(self, iters, **kwargs):
    '''optimize logits'''
    self.design(iters, **kwargs)

  def design_soft(self, iters, **kwargs):
    ''' optimize softmax(logits/temp)'''
    self.design(iters, soft=True, **kwargs)
  
  def design_hard(self, iters, **kwargs):
    ''' optimize argmax(logits)'''
    self.design(iters, soft=True, hard=True, **kwargs)

  def design_2stage(self, soft_iters=100, temp_iters=100, hard_iters=50,
                    temp=1.0, dropout=True, gumbel=False, **kwargs):
    '''two stage design (soft→hard)'''
    self.design(soft_iters, soft=True, temp=temp, dropout=dropout, gumbel=gumbel, **kwargs)
    self.design(temp_iters, soft=True, temp=temp, dropout=dropout, gumbel=False,  e_temp=1e-2, **kwargs)
    self.design(hard_iters, soft=True, temp=1e-2, dropout=False,   gumbel=False,  hard=True, save_best=True, **kwargs)

  def design_3stage(self, soft_iters=300, temp_iters=100, hard_iters=50, 
                    temp=1.0, dropout=True, gumbel=False, **kwargs):
    '''three stage design (logits→soft→hard)'''
    self.design(soft_iters, e_soft=True, temp=temp, dropout=dropout, gumbel=gumbel, **kwargs)
    self.design(temp_iters, soft=True,   temp=temp, dropout=dropout, gumbel=False, e_temp=1e-2,**kwargs)
    self.design(hard_iters, soft=True,   temp=1e-2, dropout=False,   gumbel=False, hard=True, save_best=True, **kwargs)  
    
####################################################
# AF_UTILS - various utils (save, plot, etc)
####################################################
class _af_utils:
  def _save_results(self, save_best=False, verbose=True):
    '''save the results and update trajectory'''

    # save best result
    if save_best and self._loss < self._best_loss:
      self._best_loss = self._loss
      self._best_outs = self._outs
    
    # compile losses
    self._losses.update({"model":self._model_num,"loss":self._loss, 
                        **{k:self.opt[k] for k in ["soft","hard","temp"]}})
    
    if self.protocol == "fixbb":
      # compute sequence recovery
      _aatype = self._outs["seq"].argmax(-1)
      L = min(_aatype.shape[-1], self._wt_aatype.shape[-1])
      self._losses["seqid"] = (_aatype[...,:L] == self._wt_aatype[...,:L]).mean()

    # save losses
    self.losses.append(self._losses)

    # print losses  
    if verbose:
      def to_str(ks, f=2):
        out = []
        for k in ks:
          if k in self._losses and ("rmsd" in k or self.opt["weights"].get(k,True)):
            out.append(f"{k}")
            if f is None: out.append(f"{self._losses[k]}")
            else: out.append(f"{self._losses[k]:.{f}f}")
        return out
      out = [to_str(["model","recycles"],None),
             to_str(["soft","temp","seqid","loss"]),
             to_str(["msa_ent","plddt","pae","helix","con",
                     "i_pae","i_con",
                     "sc_fape","sc_rmsd",
                     "dgram_cce","fape","rmsd"])]
      print_str = " ".join(sum(out,[]))
      print(f"{self._k}\t{print_str}")

    # save trajectory
    ca_xyz = self._outs["final_atom_positions"][:,1,:]
    traj = {"xyz":ca_xyz,"plddt":self._outs["plddt"],"seq":self._outs["seq_pseudo"]}
    if "pae" in self._outs: traj.update({"pae":self._outs["pae"]})
    for k,v in traj.items(): self._traj[k].append(np.array(v))

  def get_seqs(self):
    outs = self._outs if self._best_outs is None else self._best_outs
    outs = jax.tree_map(lambda x:np.asarray(x), outs)
    x = np.array(outs["seq"]).argmax(-1)
    return ["".join([order_restype[a] for a in s]) for s in x]
  
  def get_loss(self, x="loss"):
    '''output the loss (for entire trajectory)'''
    return np.array([float(loss[x]) for loss in self.losses])

  def save_pdb(self, filename=None):
    '''save pdb coordinates'''
    outs = self._outs if self._best_outs is None else self._best_outs
    outs = jax.tree_map(lambda x:np.asarray(x), outs)
    aatype = outs["seq"].argmax(-1)[0]
    if self.protocol == "binder":
      aatype_target = self._batch["aatype"][:self._target_len]
      aatype = np.concatenate([aatype_target,aatype])
    if self.protocol in ["fixbb","hallucination"] and self._copies > 1:
      aatype = np.concatenate([aatype] * self._copies)
    p = {"residue_index":self._inputs["residue_index"][0],
          "aatype":aatype,
          "atom_positions":outs["final_atom_positions"],
          "atom_mask":outs["final_atom_mask"]}
    b_factors = outs["plddt"][:,None] * p["atom_mask"]
    p = protein.Protein(**p,b_factors=b_factors)
    pdb_lines = protein.to_pdb(p)
    if filename is None:
      return pdb_lines
    else:
      with open(filename, 'w') as f: f.write(pdb_lines)
  
  #-------------------------------------
  # plotting functions
  #-------------------------------------
  def animate(self, s=0, e=None, dpi=100):
    sub_traj = {k:v[s:e] for k,v in self._traj.items()}
    if self.protocol == "fixbb":
      pos_ref = self._batch["all_atom_positions"][:,1,:]
      length = self._len if self._copies > 1 else None
      return make_animation(**sub_traj, pos_ref=pos_ref, length=length, dpi=dpi)
    
    if self.protocol == "binder":
      outs = self._outs if self._best_outs is None else self._best_outs
      pos_ref = outs["final_atom_positions"][:,1,:]
      return make_animation(**sub_traj, pos_ref=pos_ref,
                            length=self._target_len, dpi=dpi)

    if self.protocol in ["hallucination","partial"]:
      outs = self._outs if self._best_outs is None else self._best_outs
      pos_ref = outs["final_atom_positions"][:,1,:]
      length = self._len if self._copies > 1 else None
      return make_animation(**sub_traj, pos_ref=pos_ref, length=length, dpi=dpi)

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
  
  def plot_traj(self, dpi=100):
    fig = plt.figure(figsize=(5,5), dpi=dpi)
    gs = GridSpec(4,1, figure=fig)
    ax1 = fig.add_subplot(gs[:3,:])
    ax2 = fig.add_subplot(gs[3:,:])
    ax1_ = ax1.twinx()
    
    if self.protocol == "fixbb":
      rmsd = self.get_loss("rmsd")
      for k in [0.5,1,2,4,8,16,32]:
        ax1.plot([0,len(rmsd)],[k,k],color="lightgrey")
      ax1.plot(rmsd,color="black")
      ax1_.plot(self.get_loss("seqid"),color="green",label="seqid")
      # axes labels
      ax1.set_yscale("log")
      ticks = [0.25,0.5,1,2,4,8,16,32,64]
      ax1.set(xticks=[])
      ax1.set_yticks(ticks);ax1.set_yticklabels(ticks)
      ax1.set_ylabel("RMSD",color="black");ax1_.set_ylabel("seqid",color="green")
      ax1.set_ylim(0.25,64)
      ax1_.set_ylim(0,0.4)
      # extras
      if "soft" in self.losses[0]:
        ax2.plot(self.get_loss("soft"),color="yellow",label="soft")
      if "temp" in self.losses[0]:
        ax2.plot(self.get_loss("temp"),color="orange",label="temp")
      if "hard" in self.losses[0]:
        ax2.plot(self.get_loss("hard"),color="red",label="hard")
      ax2.set_ylim(-0.1,1.1)
      ax2.set_xlabel("iterations")
      ax2.legend(loc='center left')
    else:
      print("TODO")
    plt.show()
    
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

################################################################
# MK_DESIGN_MODEL - initialize model, and put it all together
################################################################

class mk_design_model(_af_init, _af_loss, _af_design, _af_utils):
  def __init__(self, protocol="fixbb", num_seq=1,
               num_models=1, model_mode="sample", model_parallel=False,
               num_recycles=0, recycle_mode="average",
               use_templates=False,
               use_template_tor=True):

    # decide if templates should be used
    if protocol == "binder": use_templates = True

    self.protocol = protocol
    
    self.args = {"num_seq":num_seq, "use_templates":use_templates,
                 "model_mode":model_mode, "model_parallel": model_parallel,
                 "recycle_mode":recycle_mode}
    
    self._default_opt = {"temp":1.0, "soft":0.0, "hard":0.0,"gumbel":False,
                         "dropout":True, "dropout_scale":1.0,
                         "recycles":num_recycles, "models":num_models,
                         "con":  {"num":2, "cutoff":14.0, "binary":False, "seqsep":9},
                         "i_con":{"num":1, "cutoff":20.0, "binary":False},
                         "bias":np.zeros(20),
                         "template_aatype":21, "template_dropout":0.0}

    self._default_weights = {"msa_ent":0.0, "helix":0.0,
                             "plddt":0.01, "pae":0.01}

    # setup which model params to use
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

    # number of sequences
    if use_templates:
      cfg.data.eval.max_templates = 1
      cfg.data.eval.max_msa_clusters = num_seq + 1
      cfg.model.embeddings_and_evoformer.template.embed_torsion_angles = use_template_tor

    else:
      cfg.data.eval.max_msa_clusters = num_seq

    cfg.data.common.max_extra_msa = 1
    cfg.data.eval.masked_msa_replace_fraction = 0

    # number of recycles
    if recycle_mode == "average":
      cfg.model.num_recycle = 0
      cfg.data.common.num_recycle = 0
    else:
      cfg.model.num_recycle = num_recycles
      cfg.data.common.num_recycle = num_recycles

    # backprop through recycles
    if recycle_mode == "add_prev":
      cfg.model.add_prev = True
    if recycle_mode == "backprop":
      cfg.model.backprop_recycle = True
      cfg.model.embeddings_and_evoformer.backprop_dgram = True
      cfg.model.embeddings_and_evoformer.backprop_dgram_temp = 2.0

    self._config = cfg

    # setup model
    self._model_params = [data.get_model_haiku_params(model_name=model_name, data_dir=".")]
    self._runner = model.RunModel(self._config, self._model_params[0], is_training=True)

    # load the other model_params
    if use_templates:
      model_names = ["model_2_ptm"]
    else:
      model_names = ["model_1_ptm","model_2_ptm","model_4_ptm","model_5_ptm"]


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
    if protocol == "fixbb":           self.prep_inputs = self._prep_fixbb
    if protocol == "hallucination":   self.prep_inputs = self._prep_hallucination
    if protocol == "binder":          self.prep_inputs = self._prep_binder
    if protocol == "partial":         self.prep_inputs = self._prep_partial
