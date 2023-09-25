import jax
import jax.numpy as jnp
import numpy as np
import re

from colabdesign.af.alphafold.common import protein, residue_constants
from colabdesign.af.alphafold.model.tf import shape_placeholders
from colabdesign.af.alphafold.model import config

from colabdesign.shared.protein import _np_get_cb, pdb_to_string
from colabdesign.shared.prep import prep_pos
from colabdesign.shared.utils import copy_dict
from colabdesign.shared.model import order_aa

resname_to_idx = residue_constants.resname_to_idx
idx_to_resname = dict((v,k) for k,v in resname_to_idx.items())

#################################################
# AF_PREP - input prep functions
#################################################
class _af_prep:

  def _prep_model(self, **kwargs):
    '''prep model'''
    if not hasattr(self,"_model") or self._cfg != self._model["runner"].config:
      self._cfg.model.global_config.subbatch_size = None
      self._model = self._get_model(self._cfg)
      if sum(self._lengths) > 384:
        self._cfg.model.global_config.subbatch_size = 4
        self._model["fn"] = self._get_model(self._cfg)["fn"]

    self.restart(**kwargs)

  def _prep_contigs(self,
    contigs,
    pdb_filename=None,
    copies=1,
    ignore_missing=True,
    parse_as_homooligomer=False,
    partial_loss=True,
    **kwargs):
    '''prep inputs for partial hallucination protocol'''
  
    # position constraints
    self._opt["template"]["rm_ic"] = kwargs.pop("rm_template_ic",False)
    pos_cst = {}
    keys = ["hotspot", "fix_pos", "rm_template", "rm_template_seq", "rm_template_sc"]
    for k in keys:
      v = kwargs.pop(k, None)
      if v is not None:
        if isinstance(v,str):
          # to be parsed later with contigs
          pos_cst[k] = [v,[]]

        elif isinstance(v,bool):
          # set all positions to True/False
          if "template" in k:
            k_ = {"rm_template":"rm","rm_template_seq":"rm_seq","rm_template_sc":"rm_sc"}[k]
            self._opt["template"][k_] = v
          elif k == "fix_pos":
            kwargs["fix_all_pos"] = True

    # parse contigs
    if isinstance(contigs,str): contigs = contigs.replace(",",":").split(":")
    if isinstance(contigs,int): contigs = [contigs]
    
    chain_info = {}
    length, batch = [],[]
    unsupervised, fix_pos_list, hotspot_list = [],[],[]
    total_len = 0

    # for each contig
    for contig in contigs:
      if isinstance(contig,int): contig = str(contig)
      contig_len = 0
      contig_batch = []

      # for each segment in contig
      for seg in contig.split("/"):
        if seg[0].isalpha():
          # supervised
          chain = seg[0]
          if chain not in chain_info:
            chain_info[chain] = prep_pdb(pdb_filename, chain=chain, ignore_missing=ignore_missing)
          
          pos = prep_pos(seg, **chain_info[chain]["idx"])["pos"]
          seg_len = len(pos)
          contig_batch.append(jax.tree_map(lambda x:x[pos], chain_info[chain]["batch"]))
          unsup = (contig_batch[-1]["aatype"] == -1).tolist()
          unsupervised += unsup
          
          # parse position constraints
          pos_list = pos.tolist()
          for k,(x,_) in pos_cst.items():
            x_pos_list = prep_pos(x, **chain_info[chain]["idx"], error_chk=False)["pos"].tolist()
            for p in x_pos_list:
              if p in pos_list:
                i = pos_list.index(p)
                if k != "fix_pos" or not unsup[i]:
                  pos_cst[k][1].append(i + total_len)
        else:
          # unsupervised
          if "-" in seg: seg = np.random.randint(*(int(x) for x in seg.split("-")))
          seg_len = int(seg)
          contig_batch.append({'aatype': np.full(seg_len,-1),
                               'all_atom_positions': np.zeros((seg_len,37,3)),
                               'all_atom_mask': np.zeros((seg_len,37)),
                               'residue_index': np.arange(seg_len)})
          unsupervised += [True] * seg_len
        contig_len += seg_len
        total_len += seg_len
      
      # adjust residue index to be continious
      for b in range(len(contig_batch)):
        offset = -contig_batch[b]["residue_index"][0]
        if b > 0:
          offset += contig_batch[b-1]["residue_index"][-1]
        contig_batch[b]["residue_index"] += offset + 1
      
      length.append(contig_len)
      batch.append(contig_batch)

    # cst
    pos_cst = {k:np.array(v) for k,(x,v) in pos_cst.items() if len(v) > 0}

    # concatenate batch
    batch = jax.tree_map(lambda *x:np.concatenate(x),*sum(batch,[]))
    
    # propagate batch over copies
    self._args["copies"] = copies
    if copies > 1:
      if kwargs.pop("homooligomer",parse_as_homooligomer):
        assert len(length) % copies == 0
        units = len(length) // copies
        length = length[:units]
        L = sum(length)
        interchain_mask = np.ones((copies*L,copies*L))
        pos_cst = jax.tree_map(lambda x:x[x<L], pos_cst)

      else:
        unsupervised = unsupervised * copies
        batch = [batch]*copies
        batch = jax.tree_map(lambda *x:np.concatenate(x),*batch)
        L = sum(length)
        interchain_mask = np.full((copies,copies,L,L),False)
        interchain_mask[np.diag_indices(copies)] = True
        interchain_mask = interchain_mask.swapaxes(1,2).reshape(copies*L,copies*L)

      residue_index = batch["residue_index"][:sum(length)]
    else:
      L = sum(length)
      interchain_mask = np.full((L,L),True)
      residue_index = batch["residue_index"]

    # encode chains
    self._len = sum(length)
    self._lengths, chain_enc = get_chain_enc(copies, length, residue_index)
    self._wt_aatype = batch["aatype"][:self._len]
    num_res = sum(self._lengths)
        
    # configure input features
    self._inputs = prep_input_features(L=num_res, N=self._num, C=self._args["copies"],
      T=self._args["num_templates"], one_hot_msa=not self._args["optimize_seq"])
    self._inputs.update(chain_enc)
    self._inputs["wt_aatype"] = batch["aatype"]
    self._inputs["batch"] = batch
    
    # define masks
    unsupervised = np.array(unsupervised)
    unsupervised_2d = np.logical_or(unsupervised[:,None], unsupervised[None,:])

    fix_pos = np.full(num_res,kwargs.pop("fix_all_pos",False))
    if "fix_pos" in pos_cst:
      fix_pos[pos_cst.pop("fix_pos")] = True

    self._inputs.update({
      "unsupervised":    unsupervised,
      "unsupervised_2d": unsupervised_2d,
      "supervised":      unsupervised == False,
      "supervised_2d":   unsupervised_2d == False,
      "interchain_mask": interchain_mask,
      "fix_pos":         fix_pos
    })
    
    # set positional constraints
    for k,p in pos_cst.items():
      m = np.full(num_res,False)
      m[p] = True
      self._inputs[k] = m

    # decide on weights
    self.set_opt(weights=0, partial_loss=partial_loss, set_defaults=True)
    unsup_mask = np.logical_or(unsupervised_2d, interchain_mask == False)
    intra_mask = self._inputs["asym_id"][:,None] == self._inputs["asym_id"][None,:]
    self.set_weights(
      con=np.logical_and(unsup_mask,intra_mask).any(),
      i_con=np.logical_and(unsup_mask,intra_mask == False).any(),
      dgram_cce=self._inputs["supervised_2d"].any(), set_defaults=True)

    # prep model
    self._prep_model(**kwargs)

  def _prep_hallucination(self, length=100, copies=1, **kwargs):
    # prep model
    if isinstance(length,int): length = [length]
    self._prep_contigs([str(L) for L in length], copies=copies, **kwargs)

  def _prep_fixbb(self, 
    pdb_filename, 
    chain="A",
    copies=1,
    fix_pos=None,
    parse_as_homooligomer=False,
    ignore_missing=True,
    partial_loss=False,
    **kwargs):
    
    # parsing inputs
    chains, contigs = parse_chain(chain)

    # prep model
    self._prep_contigs(contigs, pdb_filename,
      fix_pos=fix_pos, ignore_missing=ignore_missing, 
      copies=copies, parse_as_homooligomer=parse_as_homooligomer, 
      partial_loss=partial_loss, **kwargs)

  def _prep_partial(self,
    pdb_filename,
    chain="A",
    length=None,
    copies=1, 
    pos=None, 
    fix_pos=None,
    parse_as_homooligomer=False,
    ignore_missing=True,
    partial_loss=True,
    **kwargs):

    # parsing inputs
    chains = chain.split(",") if isinstance(chain,str) else chain
    
    # prep contigs
    contigs = [[] for c in chains]
    if pos is not None:
      for p in pos.split(","):
        if p[0].isalpha():
          if p[0] in chain:
            contigs[chains.index(p[0])].append(p)
        else:
          contigs[0].append(chains[0]+p)
    for n,chain in enumerate(chains):
      if len(contigs[n]) == 0:
        contigs[n].append(chain)
      else:
        # fill in missing regions
        num_segments = len(contigs[n])
        if num_segments > 1:
          new_contig = [contigs[n][0]]
          for i in range(num_segments-1):
            (a, b) = (contigs[n][i], contigs[n][i+1])
            if a[0] == b[0] and len(a) > 1 and len(b) > 1:
              (a_end, b_start) = (int(a.split("-")[-1]), int(b.split("-")[0]))
              gap = (b_start - a_end) - 1
              if gap > 0:
                new_contig.append(gap)
            new_contig.append(b)
          contigs[n] = new_contig

    contigs = ["/".join(contig) for contig in contigs]

    print(f"WARNING: 'partial' protocol is being deprecated in favor of a unified 'contigs' protocol!")
    print(f"You can use the 'contigs' protocol to get the same result (but now with more flexible control):")
    print(f"model = mk_af_model('contigs')")
    print(f"model.prep_inputs(contigs='{':'.join(contigs)}', pdb_filename='{pdb_filename}')")
    
    # prep model
    self._prep_contigs(contigs, pdb_filename,
      fix_pos=fix_pos, ignore_missing=ignore_missing, 
      copies=copies, parse_as_homooligomer=parse_as_homooligomer, 
      partial_loss=partial_loss, **kwargs)

  def _prep_binder(self,
    pdb_filename,
    target_chain="A",
    binder_len=50,
    binder_chain=None,
    hotspot=None,
    fix_pos=None,
    target_flexible=False,
    ignore_missing=True,
    **kwargs):

    # parsing inputs
    target_chains, target_contigs = parse_chain(target_chain)
    binder_chains, binder_contigs = parse_chain(binder_chain)
    assert len(binder_contigs) > 0 or isinstance(binder_len,int)
    if len(binder_contigs) == 0:
      if isinstance(binder_len,int): binder_len = [binder_len]
      for i in binder_len:
        binder_contigs.append(str(i))

    # define contigs and fix_pos
    contigs_list = target_contigs + binder_contigs
    fix_pos_list = target_chains    
    if len(binder_chains) == 0:
      redesign = False
    else:
      redesign = True
      if fix_pos is not None:
        if isinstance(fix_pos,str):
          fix_pos = fix_pos.split(",")
        for a in fix_pos:
          if a[0].isalpha():
            if a[0].isalpha() in binder_chains:
              fix_pos_list.append(a)
          else:
            fix_pos_list.append(binder_chains[0]+a)
    
    contigs = ":".join(contigs_list)
    fix_pos = ",".join(fix_pos_list)

    # define hotspot
    hotspot_list = []
    if hotspot is not None:
      if isinstance(hotspot,str):
        hotspot = hotspot.split(",")
      for a in hotspot:
        if a[0].isalpha():
          if a[0].isalpha() in target_chains:
            hotspot_list.append(a)
        else:
          hotspot_list.append(target_chains[0]+a)

    if len(hotspot_list) > 0:
      hotspot = ",".join(hotspot_list)

    template_args = {
      "rm_target":False, "rm_target_seq":False, "rm_target_sc":False,
      "rm_binder":True,  "rm_binder_seq":True,  "rm_binder_sc":True,
      "rm_template_ic":True
    }
    template_args.update({k:kwargs.pop(k,v) for k,v in template_args.items()})
    if target_flexible: template_args.update({"rm_target_seq":True, "rm_target_sc":True})

    self._prep_contigs(contigs,
      pdb_filename=pdb_filename,
      fix_pos=fix_pos,
      hotspot=hotspot,
      ignore_missing=ignore_missing,
      **kwargs)

    # adjust weights
    if redesign:
      self.set_weights(i_con=0, con=0, plddt=0, dgram_cce=1)
    else:
      self.set_weights(i_con=1, con=0, plddt=0.1, dgram_cce=0)

    # adjust masks
    num_res = sum(self._lengths)
    self._target_len = tL = sum(self._lengths[:len(target_contigs)])
    self._binder_len = bL = num_res - tL

    if hotspot is None:
      self._inputs["hotspot"] = np.array([0]*tL+[1]*bL)
    self._inputs["supervised"][:tL] = 0
    self._inputs["supervised_2d"][:tL,:tL] = 0

    if template_args["rm_template_ic"]:
      self._inputs["interchain_mask"][:tL,tL:] = False
      self._inputs["interchain_mask"][tL:,:tL] = False
    for k in ["","_seq","_sc"]:
      m = np.full(num_res,0)
      if template_args[f"rm_target{k}"]: m[:tL] = True
      if template_args[f"rm_binder{k}"]: m[tL:] = True
      self._inputs[f"rm_template{k}"] = m

#######################
# utils
#######################

def parse_chain(chain):
  if chain is None:
    return [],[]
  if isinstance(chain,str):
    chain = chain.replace(",",":").split(":")
  chains, contigs = [],[]
  for c in chain:
    contigs.append(c)
    for seg in c.split("/"):
      if seg[0].isalpha():
        if seg[0] not in chains:
          chains.append(seg[0])
  return chains, contigs

def prep_pdb(pdb_filename, chain=None,
             offsets=None, lengths=None,
             ignore_missing=False,
             offset_index=False, auth_chains=True):
  '''extract features from pdb'''

  def add_cb(batch):
    '''add missing CB atoms based on N,CA,C'''
    p,m = batch["all_atom_positions"], batch["all_atom_mask"]
    atom_idx = residue_constants.atom_order
    atoms = {k:p[...,atom_idx[k],:] for k in ["N","CA","C"]}
    cb = atom_idx["CB"]
    cb_atoms = _np_get_cb(**atoms, use_jax=False)
    cb_mask = np.prod([m[...,atom_idx[k]] for k in ["N","CA","C"]],0)
    batch["all_atom_positions"][...,cb,:] = np.where(m[:,cb,None], p[:,cb,:], cb_atoms)
    batch["all_atom_mask"][...,cb] = (m[:,cb] + cb_mask) > 0
    return {"atoms":batch["all_atom_positions"][:,cb],"mask":cb_mask}

  chains = chain.split(",") if isinstance(chain,str) else chain
  o,last = [],0
  residue_idx, chain_idx = [],[]
  full_lengths = []

  # go through each defined chain  
  for n,chain in enumerate(chains):
    pdb_str = pdb_to_string(pdb_filename, chains=chain, models=[1], auth_chains=auth_chains)
    protein_obj = protein.from_pdb_string(pdb_str) #, chain_id=chain)
    batch = {'aatype': protein_obj.aatype,
             'all_atom_positions': protein_obj.atom_positions,
             'all_atom_mask': protein_obj.atom_mask,
             'residue_index': protein_obj.residue_index}

    cb_feat = add_cb(batch) # add in missing cb (in the case of glycine)
    
    im = ignore_missing[n] if isinstance(ignore_missing,list) else ignore_missing
    if im:
      r = batch["all_atom_mask"][:,residue_constants.atom_order["CA"]] == 1
      batch = jax.tree_map(lambda x:x[r], batch)
      residue_index = batch["residue_index"] + last

    else:
      # pad values
      offset = 0 if offsets is None else (offsets[n] if isinstance(offsets,list) else offsets)
      r = offset + (protein_obj.residue_index - protein_obj.residue_index.min())
      length = (r.max()+1) if lengths is None else (lengths[n] if isinstance(lengths,list) else lengths)    
      def scatter(x, value=0):
        shape = (length,) + x.shape[1:]
        y = np.full(shape, value, dtype=x.dtype)
        y[r] = x
        return y

      batch.update({"aatype":scatter(batch["aatype"],-1),
        "all_atom_positions":scatter(batch["all_atom_positions"]),
        "all_atom_mask":scatter(batch["all_atom_mask"])})
      batch["residue_index"] = np.arange(batch["residue_index"][0], batch["residue_index"][-1] + 1)
      
      residue_index = np.arange(length) + last
    
    if offset_index:
      last = residue_index[-1] + 50
    else:
      last = 0
    o.append({"batch":batch,
              "residue_index": residue_index,
              "cb_feat":cb_feat})
    
    residue_idx.append(batch["residue_index"])
    chain_idx.append([chain] * len(residue_idx[-1]))
    full_lengths.append(len(residue_index))

  # concatenate chains
  o = jax.tree_util.tree_map(lambda *x:np.concatenate(x,0),*o)
  
  # save original residue and chain index
  o["idx"] = {"residue":np.concatenate(residue_idx), "chain":np.concatenate(chain_idx)}
  o["lengths"] = full_lengths
  return o

def make_fixed_size(feat, num_res, num_seq=1, num_templates=1, skip_batch=False):
  '''pad input features'''
  def pad_fn(arr, pad_values, pad_value=0):
    pad_values = np.array(pad_values)
    # Handle negative padding (removal of elements)
    for idx, (pad_start, pad_end) in enumerate(pad_values):
      if len(arr.shape) > idx:  # If the dimension exists
        if pad_start < 0:  # If start padding is negative, remove elements from the start
          arr = arr[abs(pad_start):]
        if pad_end < 0:  # If end padding is negative, remove elements from the end
          arr = arr[:arr.shape[0] + pad_end]
    # Handle positive padding (addition of elements)
    pad_values = np.where(pad_values < 0, 0, pad_values)
    arr = np.pad(arr, pad_values, mode='constant', constant_values=pad_value)
    return arr

  shape_schema = {k:v for k,v in config.CONFIG.data.eval.feat.items()}
  pad_size_map = {
      shape_placeholders.NUM_RES: num_res,
      shape_placeholders.NUM_MSA_SEQ: num_seq,
      shape_placeholders.NUM_EXTRA_SEQ: 1,
      shape_placeholders.NUM_TEMPLATES: num_templates
  }  
  for k,v in feat.items():
    if k == "batch" and not skip_batch:
      feat[k] = make_fixed_size(v, num_res, num_seq, num_templates)
    elif k in shape_schema:
      shape = list(v.shape)
      schema = shape_schema[k][:len(shape)]
      pad_size = [pad_size_map.get(s2, None) or s1 for (s1, s2) in zip(shape, schema)]
      padding = [(0, p - v.shape[i]) for i, p in enumerate(pad_size)]
      feat[k] = pad_fn(v, padding)
  return feat

def prep_input_features(L, N=1, T=1, C=1, one_hot_msa=True):
  '''
  given [L]ength, [N]umber of sequences and number of [T]emplates
  return dictionary of blank features
  '''
  sub_L = L//C
  inputs = {'msa': np.zeros((N,sub_L),int) if one_hot_msa else np.zeros((N,sub_L,22)),
            'deletion_matrix': np.zeros((N,sub_L),int),
            'bias': np.zeros((L,20)),
            'seq_mask': np.full(L,True),
            'residue_index': np.arange(L),
            'asym_id': np.zeros(L,int),
            'sym_id': np.zeros(L,int),
            'entity_id': np.zeros(L,int),

            # for template inputs
            'template_aatype': np.zeros((T,L),int),
            'template_all_atom_mask': np.zeros((T,L,37)),
            'template_all_atom_positions': np.zeros((T,L,37,3)),
            'template_mask': np.full(T,True)}
  return inputs

def get_chain_enc(copies, lengths, residue_idx=None):
  if residue_idx is None:
    residue_idx = np.concatenate([np.arange(L) for L in lengths],0)
  chain_enc = {"asym_id":[],"sym_id":[],"entity_id":[],"residue_index":[]}
  asym_id = 0
  chain_lengths = []
  for sym_id in range(copies):
    for entity_id, L in enumerate(lengths):
      chain_enc["asym_id"].append(np.full(L,asym_id))
      chain_enc["sym_id"].append(np.full(L,sym_id))
      chain_enc["entity_id"].append(np.full(L,entity_id))
      chain_lengths.append(L)
      asym_id += 1
    chain_enc["residue_index"].append(residue_idx)
  for k in chain_enc.keys():
    chain_enc[k] = np.concatenate(chain_enc[k])
  return chain_lengths, chain_enc