from af_backprop.utils import *
from af_backprop.utils import _np_get_6D_loss, _np_rmsd

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

  def _get_fixbb_loss(self, outputs, opt):
    losses = {"fape":      get_fape_loss(self._batch, outputs, model_config=self._config),
              "dgram_cce": get_dgram_loss(self._batch, outputs, model_config=self._config, copies=self._copies),
              "rmsd":      get_rmsd_loss_w(self._batch, outputs, copies=self._copies)}
    return {"losses":losses, "aux":{}}

  def _get_partial_loss(self, outputs, opt):
    '''compute loss for partial hallucination protocol'''
    losses, aux = {},{}
    pos = opt["pos"]
    _config = self._config.model.heads.structure_module

    # dgram
    dgram = outputs["distogram"]["logits"][:,pos][pos,:]
    losses["dgram_cce"] = get_dgram_loss(self._batch, outputs, model_config=self._config, logits=dgram)

    # 6D loss
    true = self._batch["all_atom_positions"]
    pred = outputs["structure_module"]["final_atom_positions"][pos]
    losses["6D"] = _np_get_6D_loss(true, pred, use_theta=not self.args["sidechain"])

    # rmsd
    losses["rmsd"] = _np_rmsd(true[:,1], pred[:,1])
    
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
    aux.update({"dgram":dgram,
                "dgram_bins":dgram_bins})

    # contact
    c,ic = opt["con"],opt["i_con"]
    def get_con(x, cutoff, binary=True, entropy=True):
      bins = dgram_bins < cutoff
      
      px = jax.nn.softmax(x)
      
      con_loss_bin = 1 - (bins * px).sum(-1)
      con_loss_cat = 1 - (bins * px).max(-1)
      con_loss = jnp.where(binary, con_loss_bin, con_loss_cat)
      
      # binary/cateogorical cross-entropy
      con_loss_bin_ent = -jnp.log((bins * px + 1e-8).sum(-1))      
      con_loss_cat_ent = -(jax.nn.softmax(x - 1e7 * (1-bins)) * jax.nn.log_softmax(x)).sum(-1)
      con_loss_ent = jnp.where(binary, con_loss_bin_ent, con_loss_cat_ent)
      
      return jnp.where(entropy, con_loss_ent, con_loss)
    
    def get_con_score(dgram,c):
      def set_diag(x, k, val=0.0):
        L = x.shape[-1]
        mask = jnp.abs(jnp.arange(L)[:,None] - jnp.arange(L)[None,:]) < k
        return jnp.where(mask, val, x)      
      
      def min_k(x, k=1):
        x = jnp.sort(x)
        mask = jnp.arange(x.shape[-1]) < k
        return jnp.where(mask, x, 0.0).sum(-1)/mask.sum(-1)
      
      if "seqsep" in c:
        x = get_con(dgram,c["cutoff"],c["binary"],c["entropy"])
        x_diag_ent = set_diag(x,c["seqsep"],1e8)
        x_diag = set_diag(x,c["seqsep"],1.0)
        x = jnp.where(c["entropy"],x_diag_ent,x_diag)
      else:
        x = get_con(dgram,c["cutoff"],c["binary"],c["entropy"])
      return min_k(x,c["num"])
    
    def get_helix_score(dgram,c):
      x = get_con(dgram, 6.0, binary=True, entropy=c["entropy"])
      return jnp.diagonal(x,3)

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
          losses["helix"] = get_helix_score(x,c).mean()
          x = get_con_score(x,c)
          ix = get_con_score(ix,ic)

        losses.update({f"i_{k}":ix.mean(), k:x.mean()})
    
    else:
      losses.update({"con":get_con_score(dgram,c).mean(),
                     "helix":get_helix_score(dgram,c).mean(),
                     "pae":pae.mean()})
    
    return {"losses":losses, "aux":aux}
