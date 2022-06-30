import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.shared.utils import Key
from colabdesign.shared.protein import _np_get_6D_loss, _np_rmsd, _np_kabsch

from colabdesign.af.misc import update_seq, update_aatype 
from colabdesign.af.misc import get_plddt, get_fape_loss, get_dgram_loss, get_rmsd_loss_w, get_sc_rmsd
from colabdesign.af.alphafold.model import model, folding, all_atom

####################################################
# AF_LOSS - setup loss function
####################################################
class _af_loss:
  def _get_model(self):

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
      inputs["dropout_rate"] = jnp.array([opt["dropout"]]).astype(float)
      
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
                  "plddt":get_plddt(outputs),
                  "prev":outputs["prev"]})

      if self._args["debug"]:
        aux.update({"outputs":outputs, "inputs":inputs})

      #######################################################################
      # LOSS
      #######################################################################
      aux["losses"] = {}
      self._get_loss(inputs=inputs, outputs=outputs, opt=opt, aux=aux)

      if self._loss_callback is not None:
        aux["losses"].update(self._loss_callback(inputs, outputs))

      w = opt["weights"]
      loss = sum([v*w[k] if k in w else v for k,v in aux["losses"].items()])
            
      return loss, aux
    
    return jax.value_and_grad(_model, has_aux=True, argnums=0), _model

  ##################################################
  # INPUT FUNCTIONS
  ##################################################
  def _get_seq(self, params, opt, aux, key):
    '''get sequence features'''
    seq = {"input":params["seq"]}
    seq_shape = params["seq"].shape

    # shuffle msa (randomly pick which sequence is query)
    if self._num > 1:
      n = jax.random.randint(key,[],0,seq_shape[0])
      seq["input"] = seq["input"].at[0].set(seq["input"][n]).at[n].set(seq["input"][0])

    # straight-through/reparameterization
    seq["logits"] = seq["input"] * opt["alpha"] + opt["bias"]
    seq["soft"] = jax.nn.softmax(seq["logits"] / opt["temp"])
    seq["hard"] = jax.nn.one_hot(seq["soft"].argmax(-1), 20)
    seq["hard"] = jax.lax.stop_gradient(seq["hard"] - seq["soft"]) + seq["soft"]

    # create pseudo sequence
    seq["pseudo"] = opt["soft"] * seq["soft"] + (1-opt["soft"]) * seq["input"]
    seq["pseudo"] = opt["hard"] * seq["hard"] + (1-opt["hard"]) * seq["pseudo"]
    
    # fix sequence for partial hallucination
    if self.protocol == "partial":
      seq_ref = jax.nn.one_hot(self._wt_aatype,20)
      fix_seq = lambda x:jnp.where(opt["fix_seq"],x.at[:,opt["pos"],:].set(seq_ref),x)
      seq = jax.tree_map(fix_seq,seq)

    aux.update({"seq":seq, "seq_pseudo":seq["pseudo"]})
    
    # protocol specific modifications to seq features
    if self.protocol == "binder":
      # concatenate target and binder sequence
      seq_target = jax.nn.one_hot(self._batch["aatype"][:self._target_len],20)
      seq_target = jnp.broadcast_to(seq_target,(self._num, *seq_target.shape))
      seq = jax.tree_map(lambda x:jnp.concatenate([seq_target,x],1), seq)
      
    if self.protocol in ["fixbb","hallucination"] and self._copies > 1:
      seq = jax.tree_map(lambda x:expand_copies(x, self._copies, self._args["block_diag"]), seq)

    return seq

  def _update_template(self, inputs, opt, key):
    ''''dynamically update template features'''
    if self.protocol in ["partial","fixbb","binder"]:      
      L = self._batch["aatype"].shape[0]
      
      if self.protocol in ["partial","fixbb"]:
        aatype = jnp.zeros(L)
        template_aatype = jnp.broadcast_to(opt["template"]["aatype"],(L,))
      
      if self.protocol == "binder":
        if self._redesign:
          aatype = jnp.asarray(self._batch["aatype"])
          aatype = aatype.at[self._target_len:].set(0)
          template_aatype = aatype.at[self._target_len:].set(opt["template"]["aatype"])
        else:
          aatype = template_aatype = self._batch["aatype"]
      
      # get pseudo-carbon-beta coordinates (carbon-alpha for glycine)
      pb, pb_mask = model.modules.pseudo_beta_fn(aatype,
                                                 self._batch["all_atom_positions"],
                                                 self._batch["all_atom_mask"])
      
      # define template features
      template_feats = {"template_aatype": template_aatype,
                        "template_all_atom_positions": self._batch["all_atom_positions"],
                        "template_all_atom_masks": self._batch["all_atom_mask"],
                        "template_pseudo_beta": pb,
                        "template_pseudo_beta_mask": pb_mask}
      
      # protocol specific template injection
      if self.protocol == "fixbb":
        for k,v in template_feats.items():
          inputs[k] = inputs[k].at[:,0].set(v)            
          if k == "template_all_atom_masks":
            inputs[k] = inputs[k].at[:,0,:,5:].set(0)
      
      if self.protocol == "binder":
        n = self._target_len
        for k,v in template_feats.items():
          inputs[k] = inputs[k].at[:,0,:n].set(v[:n])            
          inputs[k] = inputs[k].at[:,-1,n:].set(v[n:])
          if k == "template_all_atom_masks":
            inputs[k] = inputs[k].at[:,-1,n:,5:].set(0)       
      
      if self.protocol == "partial":
        p = opt["pos"]
        for k,v in template_feats.items():
          inputs[k] = inputs[k].at[:,0,p].set(v)
          # if k == "template_aatype": v = jax.nn.one_hot(v,22)
          # inputs[k] = jnp.einsum("ij,i...->j...",p,v)[None,None]
          if k == "template_all_atom_masks":
            inputs[k] = inputs[k].at[:,0,:,5:].set(0)

    # dropout template input features
    L = inputs["template_aatype"].shape[2]
    n = self._target_len if self.protocol == "binder" else 0
    pos_mask = jax.random.bernoulli(key, 1-opt["template"]["dropout"],(L,))
    inputs["template_all_atom_masks"] = inputs["template_all_atom_masks"].at[:,:,n:].multiply(pos_mask[n:,None])
    inputs["template_pseudo_beta_mask"] = inputs["template_pseudo_beta_mask"].at[:,:,n:].multiply(pos_mask[n:])

  def _loss_hallucination(self, inputs, outputs, opt, aux):
    plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
    plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)
    aux["losses"]["plddt"] = plddt_loss.mean()
    self._get_pairwise_loss(inputs, outputs, opt, aux)

  def _loss_fixbb(self, inputs, outputs, opt, aux):
    '''get losses'''
    plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
    plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)
    self._get_pairwise_loss(inputs, outputs, opt, aux)

    copies = 1 if self._args["repeat"] else self._copies
    aatype = inputs["aatype"][0]
    aux["losses"].update({"fape": get_fape_loss(self._batch, outputs, model_config=self._config),
                          "rmsd": get_rmsd_loss_w(self._batch, outputs, copies=copies),
                          "dgram_cce": get_dgram_loss(self._batch, outputs, model_config=self._config, aatype=aatype, copies=copies),
                          "plddt":plddt_loss.mean()})

  #################################################
  # LOSS FUNCTIONS
  #################################################

  def _loss_binder(self, inputs, outputs, opt, aux):
    '''get losses'''
    plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
    plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)
    aux["losses"]["plddt"] = plddt_loss[...,self._target_len:].mean()
    self._get_pairwise_loss(inputs, outputs, opt, aux, interface=True)

    if self._redesign:
      aatype = inputs["aatype"][0]
      true = self._batch["all_atom_positions"][:,1,:]
      pred = outputs["structure_module"]["final_atom_positions"][:,1,:]
      L = self._target_len
      p = true - true[:L].mean(0)
      q = pred - pred[:L].mean(0)
      p = p @ _np_kabsch(p[:L], q[:L])
      rmsd_binder = jnp.sqrt(jnp.square(p[L:]-q[L:]).sum(-1).mean(-1) + 1e-8)
      fape = get_fape_loss(self._batch, outputs, model_config=self._config)
      aux["losses"].update({"rmsd":rmsd_binder, "fape":fape})
      
      # compute cce of binder + interface
      errors = get_dgram_loss(self._batch, outputs, model_config=self._config, aatype=aatype, copies=1, return_errors=True)
      aux["losses"]["dgram_cce"] = errors["errors"][L:,:].mean()       

  def _loss_partial(self, inputs, outputs, opt, aux):
    '''get losses'''    
    plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
    plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)
    aux["losses"]["plddt"] = plddt_loss.mean()
    self._get_pairwise_loss(inputs, outputs, opt, aux)


    def sub(x, p, axis=0):
      fn = lambda y:jnp.take(y,p,axis)
      # fn = lambda y:jnp.tensordot(p,y,(-1,axis)).swapaxes(axis,0)
      return jax.tree_map(fn, x)

    pos = opt["pos"]
    aatype = inputs["aatype"][0]
    _config = self._config.model.heads.structure_module

    # dgram
    dgram = sub(sub(outputs["distogram"]["logits"],pos),pos,1)
    if aatype is not None: aatype = sub(aatype,pos,0)
    aux["losses"]["dgram_cce"] = get_dgram_loss(self._batch, outputs, model_config=self._config,
                                         logits=dgram, aatype=aatype)

    true = self._batch["all_atom_positions"]
    pred = sub(outputs["structure_module"]["final_atom_positions"],pos)

    # rmsd
    aux["losses"]["rmsd"] = _np_rmsd(true[:,1], pred[:,1])

    # fape
    fape_loss = {"loss":0.0}
    struct = outputs["structure_module"]
    traj = {"traj":sub(struct["traj"],pos,-2)}
    folding.backbone_loss(fape_loss, self._batch, traj, _config)
    aux["losses"]["fape"] = fape_loss["loss"]

    # sidechain specific losses
    if self._args["use_sidechains"]:
      # sc_fape
      pred_pos = sub(struct["final_atom14_positions"],pos)
      sc_struct = {**folding.compute_renamed_ground_truth(self._batch, pred_pos),
                   "sidechains":{k: sub(struct["sidechains"][k],pos,1) for k in ["frames","atom_pos"]}}

      aux["losses"]["sc_fape"] = folding.sidechain_loss(self._batch, sc_struct, _config)["loss"]

      # sc_rmsd
      true_aa = self._batch["aatype"]
      true_pos = all_atom.atom37_to_atom14(self._batch["all_atom_positions"],self._batch)
      aux["losses"]["sc_rmsd"] = get_sc_rmsd(true_pos, pred_pos, true_aa)    

  def _get_pairwise_loss(self, inputs, outputs, opt, aux, interface=False):
    '''get pairwise loss features'''

    # decide on what offset to use
    if "offset" in inputs:
      offset = inputs["offset"][0]
    else:
      idx = inputs["residue_index"][0]
      offset = idx[:,None] - idx[None,:]

    # pae loss
    pae_prob = jax.nn.softmax(outputs["predicted_aligned_error"]["logits"])
    pae = (pae_prob * jnp.arange(pae_prob.shape[-1])).mean(-1)
    aux["pae"] = pae * 31.0
    
    # define distogram
    dgram = outputs["distogram"]["logits"]
    dgram_bins = jnp.append(0,outputs["distogram"]["bin_edges"])
    if not interface:
      aux["losses"].update({"con":get_con_loss(dgram, dgram_bins, offset=offset, **opt["con"]).mean(),
                            "helix":get_helix_loss(dgram, dgram_bins, offset=offset, **opt["con"]),
                            "pae":pae.mean()})
    else:
      # split pae/con into inter/intra      
      L = self._target_len if self.protocol == "binder" else self._len
      H = self._hotspot if hasattr(self,"_hotspot") else None
      
      def split_feats(v):
        '''split pairwise features into intra and inter features'''
        if v is None: return None,None
        aa,bb = v[:L,:L], v[L:,L:]
        ab = v[:L,L:] if H is None else v[H,L:]
        ba = v[L:,:L] if H is None else v[L:,H]
        abba = (ab + ba.swapaxes(0,1)) / 2
        if self.protocol == "binder":
          return bb,abba.swapaxes(0,1)
        else:
          return aa,abba

      x_offset, ix_offset = split_feats(jnp.abs(offset))
      for k,v in zip(["pae","con"], [pae,dgram]):
        x, ix = split_feats(v)        
        if k == "con":
          aux["losses"]["helix"] = get_helix_loss(x, dgram_bins, x_offset)
          # dgram, dgram_bins, cutoff=None, num=1, seqsep=0, binary=True, offset=None
          x = get_con_loss(x, dgram_bins, offset=x_offset, **opt["con"])
          ix = get_con_loss(ix, dgram_bins, offset=ix_offset, **opt["i_con"])

        aux["losses"].update({k:x.mean(),f"i_{k}":ix.mean()})

#####################################################################################
#####################################################################################

def get_pw_con_loss(dgram, dgram_bins, cutoff=None, binary=True):
  '''dgram to contacts'''
  if cutoff is None: cutoff = dgram_bins[-1]
  bins = dgram_bins < cutoff  
  px = jax.nn.softmax(dgram)
  px_ = jax.nn.softmax(dgram - 1e7 * (1-bins))        
  # binary/cateogorical cross-entropy
  con_loss_cat_ent = -(px_ * jax.nn.log_softmax(dgram)).sum(-1)
  con_loss_bin_ent = -jnp.log((bins * px + 1e-8).sum(-1))
  return jnp.where(binary, con_loss_bin_ent, con_loss_cat_ent)

def get_con_loss(dgram, dgram_bins, cutoff=None, binary=True,
                 num=1, seqsep=0, offset=None):
  '''convert distogram into contact loss'''  
  x = get_pw_con_loss(dgram, dgram_bins, cutoff, binary)  
  a,b = x.shape
  if offset is None:
    mask = jnp.abs(jnp.arange(a)[:,None] - jnp.arange(b)[None,:]) >= seqsep
  else:
    mask = jnp.abs(offset) >= seqsep
  x = jnp.sort(jnp.where(mask,x,jnp.nan))
  k_mask = (jnp.arange(b) < num) * (jnp.isnan(x) == False)    
  return jnp.where(k_mask,x,0.0).sum(-1) / (k_mask.sum(-1) + 1e-8)

def get_helix_loss(dgram, dgram_bins, offset=None):
  '''helix bias loss'''
  x = get_pw_con_loss(dgram, dgram_bins, cutoff=6.0, binary=True)
  if offset is None:
    return jnp.diagonal(x,3).mean()
  else:
    mask = offset == 3
    return jnp.where(mask,x,0.0).sum() / (mask.sum() + 1e-8)

def get_contact_map(outputs, dist=8.0):
  '''get contact map from distogram'''
  dist_logits = outputs["distogram"]["logits"]
  dist_bins = jax.numpy.append(0,outputs["distogram"]["bin_edges"])
  dist_mtx = dist_bins[dist_logits.argmax(-1)]
  return (jax.nn.softmax(dist_logits) * (dist_bins < dist)).sum(-1)

def expand_copies(x, copies, block_diag=True):
  '''
  given msa (N,L,20) expand to (N*(1+copies),L*copies,22) if block_diag else (N,L*copies,22)
  '''
  if x.shape[-1] < 22:
    x = jnp.pad(x,[[0,0],[0,0],[0,22-x.shape[-1]]])
  x = jnp.tile(x,[1,copies,1])
  if copies > 1 and block_diag:
    L = x.shape[1]
    sub_L = L // copies
    y = x.reshape((-1,1,copies,sub_L,22))
    block_diag_mask = jnp.expand_dims(jnp.eye(copies),(0,3,4))
    seq = block_diag_mask * y
    gap_seq = (1-block_diag_mask) * jax.nn.one_hot(jnp.repeat(21,sub_L),22)  
    y = (seq + gap_seq).swapaxes(0,1).reshape(-1,L,22)
    return jnp.concatenate([x,y],0)
  else:
    return x