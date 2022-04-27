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
    self._default_opt.update({"6D_weights":{}})

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
                    fix_seq=True, sidechain=False, **kwargs):
    '''prep input for partial hallucination'''
    self.args.update({"sidechain":sidechain, "fix_seq":fix_seq})
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
    
    weights = {"dgram_cce":1.0,"con":1.0, "fape":0.0, "rmsd":0.0,"6D":0.0}
    if sidechain: weights.update({"sc_fape":0.0, "sc_rmsd":0.0})
    self._default_weights.update(weights)

    self.restart(set_defaults=True, **kwargs)
