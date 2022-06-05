import jax
import jax.numpy as jnp
import numpy as np

from af.src.misc import _np_get_6D_loss, _np_rmsd
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

      # set sequence
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
      if self.protocol == "partial" and (self.args["fix_seq"] or self.args["sidechain"]):
        p = opt["pos"]
        seq_ref = jax.nn.one_hot(self._wt_aatype,20)
        if jnp.issubdtype(p.dtype, jnp.integer):
          seq = jax.tree_map(lambda x:x.at[:,p,:].set(seq_ref), seq)
        else:
          seq_ref = p.T @ seq_ref
          w = 1 - seq_ref.sum(-1, keepdims=True)
          seq = jax.tree_map(lambda x:x*w + seq_ref, seq)

      # save for aux output
      aux = {"seq":seq}

      # entropy loss for msa
      if self.args["num_seq"] > 1:
        seq_prf = seq["hard"].mean(0)
        losses["msa_ent"] = -(seq_prf * jnp.log(seq_prf + 1e-8)).sum(-1).mean()
        
      if self.protocol == "binder":
        # concatenate target and binder sequence
        seq_target = jax.nn.one_hot(self._batch["aatype"][:self._target_len],20)
        seq_target = jnp.broadcast_to(seq_target,(self.args["num_seq"],*seq_target.shape))
        seq = jax.tree_map(lambda x:jnp.concatenate([seq_target,x],1), seq)
              
      if self.protocol in ["fixbb","hallucination"] and self._copies > 1:
        seq = jax.tree_map(lambda x:jnp.concatenate([x]*self._copies,1), seq)
            
      # update sequence
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
      if "offset" in inputs:
        offset = inputs["offset"][0]
      else:
        idx = inputs["residue_index"][0]
        offset = idx[:,None] - idx[None,:]
      
      o = self._get_pairwise_loss(outputs, opt, offset=offset)
      losses.update(o["losses"]); aux.update(o["aux"])

      # confidence losses      
      plddt_prob = jax.nn.softmax(outputs["predicted_lddt"]["logits"])
      plddt_loss = (plddt_prob * jnp.arange(plddt_prob.shape[-1])[::-1]).mean(-1)

      if self.protocol == "binder":
        losses["plddt"] = plddt_loss[...,self._target_len:].mean()
      else:
        losses["plddt"] = plddt_loss.mean()

      if self.protocol == "fixbb" or (self.protocol == "binder" and self._redesign):
        o = self._get_fixbb_loss(outputs, opt, aatype=aatype)
        losses.update(o["losses"]); aux.update(o["aux"])

      if self.protocol == "partial":
        o = self._get_partial_loss(outputs, opt, aatype=aatype)
        losses.update(o["losses"]); aux.update(o["aux"])

      # weighted loss
      loss = sum([v*opt["weights"][k] for k,v in losses.items()])

      # save aux outputs
      aux.update({"final_atom_positions":outputs["structure_module"]["final_atom_positions"],
                  "final_atom_mask":outputs["structure_module"]["final_atom_mask"],
                  "plddt":get_plddt(outputs), "losses":losses})

      return loss, (aux)
    
    return jax.value_and_grad(mod, has_aux=True, argnums=0), mod

  def _get_fixbb_loss(self, outputs, opt, aatype=None):
    copies = 1 if self._repeat else self._copies
    losses = {"fape":      get_fape_loss(self._batch, outputs, model_config=self._config),
              "dgram_cce": get_dgram_loss(self._batch, outputs, model_config=self._config, aatype=aatype, copies=copies),
              "rmsd":      get_rmsd_loss_w(self._batch, outputs, copies=copies)}
    return {"losses":losses, "aux":{}}
    
  def _get_partial_loss(self, outputs, opt, aatype=None):
    '''compute loss for partial hallucination protocol'''
    def sub(x, p, axis=0):
      if jnp.issubdtype(p.dtype, jnp.integer): return jnp.take(x,p,axis)
      else: return jnp.tensordot(p,x,(-1,axis)).swapaxes(axis,0)

    losses, aux = {},{}
    pos = opt["pos"]
    _config = self._config.model.heads.structure_module

    # dgram
    dgram = sub(sub(outputs["distogram"]["logits"],pos),pos,1)
    if aatype is not None: aatype = sub(aatype,pos,0)
    losses["dgram_cce"] = get_dgram_loss(self._batch, outputs, model_config=self._config,
                                         logits=dgram, aatype=aatype)

    # 6D loss
    true = self._batch["all_atom_positions"]
    pred = sub(outputs["structure_module"]["final_atom_positions"],pos)
    losses["6D"] = _np_get_6D_loss(true, pred, use_theta=not self.args["sidechain"])

    # rmsd
    losses["rmsd"] = _np_rmsd(true[:,1], pred[:,1])
    
    # fape
    fape_loss = {"loss":0.0}
    struct = outputs["structure_module"]
    traj = {"traj":sub(struct["traj"],pos,-2)}
    folding.backbone_loss(fape_loss, self._batch, traj, _config)
    losses["fape"] = fape_loss["loss"]

    # sidechain specific losses
    if self.args["sidechain"]:

      # sc_fape
      pred_pos = sub(struct["final_atom14_positions"],pos)
      sc_struct = {**folding.compute_renamed_ground_truth(self._batch, pred_pos),
                   "sidechains":{k: sub(struct["sidechains"][k],pos,1) for k in ["frames","atom_pos"]}}
      
      losses["sc_fape"] = folding.sidechain_loss(self._batch, sc_struct, _config)["loss"]

      # sc_rmsd
      true_aa = self._batch["aatype"]
      true_pos = all_atom.atom37_to_atom14(self._batch["all_atom_positions"],self._batch)
      losses["sc_rmsd"] = get_sc_rmsd(true_pos, pred_pos, true_aa)

    return {"losses":losses, "aux":aux}

  def _get_pairwise_loss(self, outputs, opt, offset=None):
    '''compute pairwise losses (contact and pae)'''

    losses, aux = {},{}
    # pae loss
    pae_prob = jax.nn.softmax(outputs["predicted_aligned_error"]["logits"])
    pae = (pae_prob * jnp.arange(pae_prob.shape[-1])).mean(-1)
    
    # define distogram
    dgram = outputs["distogram"]["logits"]
    dgram_bins = jnp.append(0,outputs["distogram"]["bin_edges"])
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
    
    # if more than 1 chain, split pae/con into inter/intra
    if self.protocol == "binder" or (self._copies > 1 and not self._repeat):
      aux["pae"] = get_pae(outputs)

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
      for k,v in zip(["pae","con"],[pae,dgram]):
        x, ix = split_feats(v)        
        if k == "con":
          losses["helix"] = get_helix_loss(x, c, x_offset)
          x, ix = get_con_loss(x, c, x_offset), get_con_loss(ix, ic, ix_offset)
        losses.update({k:x.mean(),f"i_{k}":ix.mean()})
    else:
      if self.protocol == "binder": aux["pae"] = get_pae(outputs)
      losses.update({"con":get_con_loss(dgram, c, offset).mean(),
                     "helix":get_helix_loss(dgram, c, offset),
                     "pae":pae.mean()})
    
    return {"losses":losses, "aux":aux}
  
  def _update_template(self, inputs, opt, key):
    ''''dynamically update template features'''
    # add template features
    t = "template_"    
    if self.protocol == "partial":
      pb, pb_mask = model.modules.pseudo_beta_fn(self._batch["aatype"],
                                                 self._batch["all_atom_positions"],
                                                 self._batch["all_atom_mask"])      
      feats = {t+"aatype": self._batch["aatype"],
               t+"all_atom_positions": self._batch["all_atom_positions"],
               t+"all_atom_masks": self._batch["all_atom_mask"],
               t+"pseudo_beta": pb,
               t+"pseudo_beta_mask": pb_mask}
            
      p = opt["pos"]      
      for k,v in feats.items():
        if jnp.issubdtype(p.dtype, jnp.integer):
          inputs[k] = inputs[k].at[:,:,p].set(v)
        else:
          if k == t+"aatype": v = jax.nn.one_hot(v,22)
          inputs[k] = jnp.einsum("ij,i...->j...",p,v)[None,None]
          
    if self.protocol == "fixbb":
      inputs[t+"aatype"] = inputs[t+"aatype"].at[:].set(opt[t+"aatype"])        
      inputs[t+"all_atom_masks"] = inputs[t+"all_atom_masks"].at[...,5:].set(0.0)
    
    # dropout template input features
    L = self._len
    s = self._target_len if self.protocol == "binder" else 0
    pos_mask = jax.random.bernoulli(key, 1-opt[t+"dropout"],(L,))
    inputs[t+"all_atom_masks"] = inputs[t+"all_atom_masks"].at[...,s:,:].multiply(pos_mask[s:,None])
    inputs[t+"pseudo_beta_mask"] = inputs[t+"pseudo_beta_mask"].at[...,s:].multiply(pos_mask[s:])
