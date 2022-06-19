import jax
import jax.numpy as jnp
import numpy as np

from af.src.misc import _np_get_6D_loss, _np_rmsd, _np_kabsch
from af.src.misc import update_seq, update_aatype 
from af.src.misc import get_plddt, get_pae, get_fape_loss, get_dgram_loss, get_rmsd_loss_w, get_sc_rmsd

from alphafold.model import model, folding, all_atom

####################################################
# AF_LOSS - setup loss function
####################################################
class _af_loss:
  def _get_fn(self):
    # setup function to get gradients
    def mod(params, model_params, inputs, key, opt):

      # initialize the loss function
      losses = {}

      # get sequence
      key, key_ = jax.random.split(key)
      seq = self._get_seq(params, opt, key_)      
      aux = {"seq":seq, "seq_pseudo":seq["pseudo"]}

      # entropy loss for msa
      if self.args["num_seq"] > 1:
        seq_prf = seq["hard"].mean(0)
        losses["msa_ent"] = -(seq_prf * jnp.log(seq_prf + 1e-8)).sum(-1).mean()
      
      # protocol specific modifications to seq features
      if self.protocol == "binder":
        # concatenate target and binder sequence
        seq_target = jax.nn.one_hot(self._batch["aatype"][:self._target_len],20)
        seq_target = jnp.broadcast_to(seq_target,(self.args["num_seq"],*seq_target.shape))
        seq = jax.tree_map(lambda x:jnp.concatenate([seq_target,x],1), seq)
        
      if self.protocol in ["fixbb","hallucination"] and self._copies > 1:
        seq = jax.tree_map(lambda x:expand_copies(x, self._copies, self.args["block_diag"]), seq)
            
      # update sequence features
      seq_pssm = seq["soft"] if self.args["use_pssm"] else None
      update_seq(seq["pseudo"], inputs, seq_pssm=seq_pssm)
      
      # update amino acid sidechain identity
      B,L = inputs["aatype"].shape[:2]
      aatype = jax.nn.one_hot(seq["pseudo"][0].argmax(-1),21)
      update_aatype(jnp.broadcast_to(aatype,(B,L,21)), inputs)
      
      # update template features
      if self.args["use_templates"]:
        key, key_ = jax.random.split(key)
        self._update_template(inputs, opt, key_)
      
      # scale dropout rate
      inputs["scale_rate"] = jnp.where(opt["dropout"],jnp.full(1,opt["dropout_scale"]),jnp.zeros(1))
      
      # get outputs
      outputs = self._get_out(inputs, model_params, opt, key)
      aux["prev"] = outputs["prev"]
      values = {"losses":losses,"aux":aux,"outputs":outputs,"opt":opt}
      
      # pairwise loss
      self._get_pairwise_loss(**values, inputs=inputs)

      # confidence losses
      if self.use_struct:
        plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
        plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)
        if self.protocol == "binder":
          if self._redesign: self._get_binder_redesign_loss(**values, aatype=aatype)
          losses["plddt"] = plddt_loss[...,self._target_len:].mean()
        else:
          losses["plddt"] = plddt_loss.mean()
          
      # get protocol specific losses
      if self.protocol == "fixbb":
        self._get_fixbb_loss(**values, aatype=aatype)                  
      if self.protocol == "partial":
        self._get_partial_loss(**values, aatype=aatype)

      # weighted loss
      loss = []
      for k,v in losses.items():
        if k in opt["weights"]: loss.append(v * opt["weights"][k])
      loss = sum(loss)
      
      aux["losses"] = losses

      # add aux outputs
      if self.use_struct:
        aux.update({"final_atom_positions":outputs["structure_module"]["final_atom_positions"],
                    "final_atom_mask":outputs["structure_module"]["final_atom_mask"],
                    "plddt":get_plddt(outputs)})
      
      if "contact_map" not in aux:
        aux["contact_map"] = get_contact_map(outputs)          
        
      if self.args["debug"]: aux.update({"outputs":outputs, "inputs":inputs})
      return loss, (aux)
    
    return jax.value_and_grad(mod, has_aux=True, argnums=0), mod

  def _get_out(self, inputs, model_params, opt, key):
    '''get outputs'''
    mode = self.args["recycle_mode"]

    def recycle(prev, key):
      inputs["prev"] = prev
      outputs = self._runner.apply(model_params, key, inputs)
      prev = outputs["prev"]
      if mode != "backprop":
        prev = jax.lax.stop_gradient(prev)
      return prev, outputs

    if mode in ["backprop","add_prev"]:
      num_iters = self.opt["recycles"] + 1
      keys = jax.random.split(key, num_iters)
      _, o = jax.lax.scan(recycle, inputs["prev"], keys)
      outputs = jax.tree_map(lambda x:x[-1], o)
      
      if mode == "add_prev":
        for k in ["distogram","predicted_lddt","predicted_aligned_error"]:
          outputs[k]["logits"] = o[k]["logits"].mean(0)
    
    elif mode in ["last","sample"]:
      def body(x):
        i,prev,key = x
        key,_key = jax.random.split(key)
        prev, _ = recycle(prev, _key)          
        return (x[0]+1, prev, key)
        
      init = (0,inputs["prev"],key)
      _, prev, key = jax.lax.while_loop(lambda x: x[0] < opt["recycles"], body, init)
      _, outputs = recycle(prev, key)
    
    else:
      outputs = self._runner.apply(model_params, key, inputs)

    return outputs

  def _get_seq(self, params, opt, key):
    '''get sequence features'''
    seq = {"input":params["seq"]}
    seq_shape = params["seq"].shape

    # shuffle msa (randomly pick which sequence is query)
    if self.args["num_seq"] > 1:
      n = jax.random.randint(key,[],0,seq_shape[0])
      seq["input"] = seq["input"].at[0].set(seq["input"][n]).at[n].set(seq["input"][0])

    # straight-through/reparameterization
    seq["logits"] = 2.0 * seq["input"] + opt["bias"] + jnp.where(opt["gumbel"], jax.random.gumbel(key,seq_shape), 0.0)
    seq["soft"] = jax.nn.softmax(seq["logits"] / opt["temp"])
    seq["hard"] = jax.nn.one_hot(seq["soft"].argmax(-1), 20)
    seq["hard"] = jax.lax.stop_gradient(seq["hard"] - seq["soft"]) + seq["soft"]

    # create pseudo sequence
    seq["pseudo"] = opt["soft"] * seq["soft"] + (1-opt["soft"]) * seq["input"]
    seq["pseudo"] = opt["hard"] * seq["hard"] + (1-opt["hard"]) * seq["pseudo"]
    
    # for partial hallucination, force sequence if sidechain constraints defined
    if self.protocol == "partial" and (self.args["fix_seq"] or self.args["use_sidechains"]):
      p = opt["pos"]
      seq_ref = jax.nn.one_hot(self._wt_aatype,20)
      if jnp.issubdtype(p.dtype, jnp.integer):
        seq = jax.tree_map(lambda x:x.at[:,p,:].set(seq_ref), seq)
      else:
        # TODO confirm this is correct
        seq_ref = p.T @ seq_ref
        w = 1 - seq_ref.sum(-1, keepdims=True)
        seq = jax.tree_map(lambda x:x*w + seq_ref, seq)
    return seq
  
  def _get_fixbb_loss(self, losses, aux, outputs, opt, aatype=None):
    '''get backbone specific losses'''
    copies = 1 if self.args["repeat"] else self._copies
    if self.use_struct:    
      losses.update({"fape": get_fape_loss(self._batch, outputs, model_config=self._config),
                     "rmsd": get_rmsd_loss_w(self._batch, outputs, copies=copies)})
    losses["dgram_cce"] = get_dgram_loss(self._batch, outputs, model_config=self._config, aatype=aatype, copies=copies)

  
  def _get_binder_redesign_loss(self, losses, aux, outputs, opt, aatype=None):
    '''get binder redesign specific losses'''
    # compute rmsd of binder relative to target
    if self.use_struct:
      true = self._batch["all_atom_positions"][:,1,:]
      pred = outputs["structure_module"]["final_atom_positions"][:,1,:]
      L = self._target_len
      p = true - true[:L].mean(0)
      q = pred - pred[:L].mean(0)
      p = p @ _np_kabsch(p[:L], q[:L])
      rmsd_binder = jnp.sqrt(jnp.square(p[L:]-q[L:]).sum(-1).mean(-1) + 1e-8)
      fape = get_fape_loss(self._batch, outputs, model_config=self._config)
      losses.update({"rmsd":rmsd_binder, "fape":fape})
    
    # compute cce of binder + interface
    errors = get_dgram_loss(self._batch, outputs, model_config=self._config, aatype=aatype, copies=1, return_errors=True)
    losses["dgram_cce"] = errors["errors"][L:,:].mean()
    
  def _get_partial_loss(self, losses, aux, outputs, opt, aatype=None):
    '''compute loss for partial hallucination protocol'''
    def sub(x, p, axis=0):
      if jnp.issubdtype(p.dtype, jnp.integer): return jnp.take(x,p,axis)
      else: return jnp.tensordot(p,x,(-1,axis)).swapaxes(axis,0)

    pos = opt["pos"]
    _config = self._config.model.heads.structure_module

    # dgram
    dgram = sub(sub(outputs["distogram"]["logits"],pos),pos,1)
    if aatype is not None: aatype = sub(aatype,pos,0)
    losses["dgram_cce"] = get_dgram_loss(self._batch, outputs, model_config=self._config,
                                         logits=dgram, aatype=aatype)

    if self.use_struct:
      # 6D loss
      true = self._batch["all_atom_positions"]
      pred = sub(outputs["structure_module"]["final_atom_positions"],pos)
      losses["6D"] = _np_get_6D_loss(true, pred, use_theta=not self.args["use_sidechains"])

      # rmsd
      losses["rmsd"] = _np_rmsd(true[:,1], pred[:,1])

      # fape
      fape_loss = {"loss":0.0}
      struct = outputs["structure_module"]
      traj = {"traj":sub(struct["traj"],pos,-2)}
      folding.backbone_loss(fape_loss, self._batch, traj, _config)
      losses["fape"] = fape_loss["loss"]

      # sidechain specific losses
      if self.args["use_sidechains"]:

        # sc_fape
        pred_pos = sub(struct["final_atom14_positions"],pos)
        sc_struct = {**folding.compute_renamed_ground_truth(self._batch, pred_pos),
                     "sidechains":{k: sub(struct["sidechains"][k],pos,1) for k in ["frames","atom_pos"]}}

        losses["sc_fape"] = folding.sidechain_loss(self._batch, sc_struct, _config)["loss"]

        # sc_rmsd
        true_aa = self._batch["aatype"]
        true_pos = all_atom.atom37_to_atom14(self._batch["all_atom_positions"],self._batch)
        losses["sc_rmsd"] = get_sc_rmsd(true_pos, pred_pos, true_aa)

  def _get_pairwise_loss(self, losses, aux, outputs, opt, inputs=None):
    '''compute pairwise losses (contact and pae)'''

    if inputs is None:
      offset = None
    else:
      if "offset" in inputs:
        offset = inputs["offset"][0] if "offset" in inputs else None
      else:
        idx = inputs["residue_index"][0]
        offset = idx[:,None] - idx[None,:]

    # pae loss
    if self.use_struct:
      pae_prob = jax.nn.softmax(outputs["predicted_aligned_error"]["logits"])
      pae = (pae_prob * jnp.arange(pae_prob.shape[-1])).mean(-1)
    
    # define distogram
    dgram = outputs["distogram"]["logits"]
    dgram_bins = jnp.append(0,outputs["distogram"]["bin_edges"])
    if self.args["debug"]:
      aux.update({"dgram":dgram, "dgram_bins":dgram_bins})

    # contact
    def get_pw_con_loss(dgram, cutoff, binary=True, entropy=True):
      '''convert distogram into pairwise contact loss'''
      bins = dgram_bins < cutoff
      
      px = jax.nn.softmax(dgram)
      px_ = jax.nn.softmax(dgram - 1e7 * (1-bins))
      
      # con_loss_cat = (bins * px_ * (1-px)).sum(-1) 
      con_loss_cat = 1 - (bins * px).max(-1)
      con_loss_bin = 1 - (bins * px).sum(-1)
      con_loss = jnp.where(binary, con_loss_bin, con_loss_cat)
      
      # binary/cateogorical cross-entropy
      con_loss_cat_ent = -(px_ * jax.nn.log_softmax(dgram)).sum(-1)
      con_loss_bin_ent = -jnp.log((bins * px + 1e-8).sum(-1))      
      con_loss_ent = jnp.where(binary, con_loss_bin_ent, con_loss_cat_ent)
      
      return jnp.where(entropy, con_loss_ent, con_loss)
    
    def get_con_loss(dgram, c, offset=None):
      '''convert distogram into contact loss'''  
      def min_k(x, k=1, seqsep=0):
        a,b = x.shape
        if offset is None:
          mask = jnp.abs(jnp.arange(a)[:,None] - jnp.arange(b)[None,:]) >= seqsep
        else:
          mask = jnp.abs(offset) >= seqsep
          
        x = jnp.sort(jnp.where(mask,x,jnp.nan))
        k_mask = (jnp.arange(b) < k) * (jnp.isnan(x) == False)
        
        return jnp.where(k_mask,x,0.0).sum(-1) / (k_mask.sum(-1) + 1e-8)

      x = get_pw_con_loss(dgram,c["cutoff"],c["binary"],c["entropy"])
      if "seqsep" in c:
        return min_k(x, c["num"], seqsep=c["seqsep"])
      else:
        return min_k(x, c["num"])
    
    def get_helix_loss(dgram,c,offset=None):
      '''helix bias loss'''
      x = get_pw_con_loss(dgram, 6.0, binary=True, entropy=c["entropy"])
      if offset is None:
        return jnp.diagonal(x,3).mean()
      else:
        mask = offset == 3
        return jnp.where(mask,x,0.0).sum() / (mask.sum() + 1e-8)

    # get contact options
    c,ic = opt["con"],opt["i_con"]
    
    if self.use_struct:
      aux["pae"] = get_pae(outputs)
      
    if self.protocol == "binder":
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

      if self.use_struct:
        ks,vs = ["pae","con"], [pae,dgram]
      else:
        ks,vs = ["con"], [dgram]
      x_offset, ix_offset = split_feats(jnp.abs(offset))
      for k,v in zip(ks,vs):
        x, ix = split_feats(v)        
        if k == "con":
          losses["helix"] = get_helix_loss(x, c, x_offset)
          x = get_con_loss(x, c, x_offset)
          ix = get_con_loss(ix, ic, ix_offset)
        losses.update({k:x.mean(),f"i_{k}":ix.mean()})
    else:
      if self.use_struct:
        if self.protocol == "binder": aux["pae"] = get_pae(outputs)
        losses.update({"con":get_con_loss(dgram, c, offset).mean(),
                       "helix":get_helix_loss(dgram, c, offset),
                       "pae":pae.mean()})
      else:
        losses.update({"con":get_con_loss(dgram, c, offset).mean(),
                       "helix":get_helix_loss(dgram, c, offset)})
              
  def _update_template(self, inputs, opt, key):
    ''''dynamically update template features'''
    if self.protocol in ["partial","fixbb","binder"]:      
      L = self._batch["aatype"].shape[0]
      
      if self.protocol in ["partial","fixbb"]:
        aatype = jnp.zeros(L)
        template_aatype = jnp.broadcast_to(opt["template_aatype"],(L,))
      
      if self.protocol == "binder":
        if self._redesign:
          aatype = jnp.asarray(self._batch["aatype"])
          aatype = aatype.at[self._target_len:].set(0)
          template_aatype = aatype.at[self._target_len:].set(opt["template_aatype"])
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
          if jnp.issubdtype(p.dtype, jnp.integer):
            inputs[k] = inputs[k].at[:,0,p].set(v)
          else:
            if k == "template_aatype": v = jax.nn.one_hot(v,22)
            inputs[k] = jnp.einsum("ij,i...->j...",p,v)[None,None]
          if k == "template_all_atom_masks":
            inputs[k] = inputs[k].at[:,0,:,5:].set(0)

    # dropout template input features
    L = inputs["template_aatype"].shape[2]
    n = self._target_len if self.protocol == "binder" else 0
    pos_mask = jax.random.bernoulli(key, 1-opt["template_dropout"],(L,))
    inputs["template_all_atom_masks"] = inputs["template_all_atom_masks"].at[:,:,n:].multiply(pos_mask[n:,None])
    inputs["template_pseudo_beta_mask"] = inputs["template_pseudo_beta_mask"].at[:,:,n:].multiply(pos_mask[n:])

def expand_copies(x, copies, block_diag=True):
  '''given msa (N,L,20) expand to 
  if block_diag:
    (N*(1+copies),L*copies,22)
  else:
    (N,L*copies,22)
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

def get_contact_map(outputs, dist=8.0):
  '''get contact map from distogram'''
  dist_logits = outputs["distogram"]["logits"]
  dist_bins = jax.numpy.append(0,outputs["distogram"]["bin_edges"])
  dist_mtx = dist_bins[dist_logits.argmax(-1)]
  return (jax.nn.softmax(dist_logits) * (dist_bins < dist)).sum(-1)
