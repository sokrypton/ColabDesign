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

    self._opt = copy_dict(self.opt)  
    self.restart(**kwargs)

  def _prep_features(self, num_res, num_seq=None, num_templates=None):
    '''process features'''
    if num_seq is None: num_seq = self._num
    if num_templates is None: num_templates = self._args["num_templates"]
    return prep_input_features(L=num_res, N=num_seq, T=num_templates,
      one_hot_msa=not self._args["optimize_seq"])

  def _prep_contigs(self, contigs, pdb_filename=None, copies=1,
    fix_pos=None, fix_all_pos=False, ignore_missing=True, 
    parse_as_homooligomer=False, **kwargs):
    '''prep inputs for partial hallucination protocol'''

    # parse contigs
    if isinstance(contigs,str): contigs = contigs.split(":")
    chain_info = {}
    length, batch = [],[]
    unsupervised, fix_pos_list = [],[]
    total_len = 0
    
    # for each contig
    for contig in contigs:
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
          unsupervised += [False] * seg_len

          # add fixed positions
          if fix_pos is not None:
            f_pos = prep_pos(fix_pos, **chain_info[chain]["idx"], error_chk=False)["pos"]
            pos_list = pos.tolist()
            for p in f_pos.tolist():
              if p in pos: fix_pos_list.append(pos_list.index(p) + total_len)
          
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

    fix_pos = np.array(fix_pos_list)

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
        if len(fix_pos_list) > 0:
          fix_pos = fix_pos[fix_pos < L]
        interchain_mask = np.ones((copies*L,copies*L))
      else:
        unsupervised = unsupervised * copies
        batch = [batch]*copies
        batch = jax.tree_map(lambda *x:np.concatenate(x),*batch)
        L = sum(length)
        interchain_mask = np.zeros((copies,copies,L,L))
        interchain_mask[np.diag_indices(copies)] = 1
        interchain_mask = interchain_mask.swapaxes(1,2).reshape(copies*L,copies*L)
      residue_index = batch["residue_index"][:sum(length)]
    else:
      L = sum(length)
      interchain_mask = np.ones((L,L))
      residue_index = batch["residue_index"]

    # encode chains
    self._len = sum(length)
    self._lengths, chain_enc = get_chain_enc(copies, length, residue_index)
    self._wt_aatype = batch["aatype"][:self._len]
        
    # configure input features
    self._inputs = self._prep_features(num_res=sum(self._lengths))
    self._inputs.update(chain_enc)
    self._inputs["batch"] = batch
    
    # define masks
    unsupervised = np.array(unsupervised)
    unsupervised_2d = np.logical_or(unsupervised[:,None],unsupervised[None,:])
    unsupervised_2d[interchain_mask == 0] = 1
    self._inputs["unsupervised"] = unsupervised
    self._inputs["unsupervised_2d"] = unsupervised_2d
    self._inputs["supervised"] = unsupervised == False
    self._inputs["supervised_2d"] = unsupervised_2d == False
    self._inputs["interchain_mask"] = interchain_mask
    
    if fix_all_pos:
      self.opt["fix_pos"] = np.arange(self._len)[unsupervised[:self_len] == False]
    elif len(fix_pos_list) > 0:
      self.opt["fix_pos"] = fix_pos

    # decide on weights
    weights = copy_dict(self.opt["weights"]) 
    losses = ["plddt","pae","i_pae","exp_res","helix","con","i_con",
      "rmsd","dgram_cce", "fape", "sc_fape","sc_chi","sc_chi_norm"]
    weights.update({k:0.0 for k in losses})
    if self._inputs["unsupervised"].sum():
      weights["con"] = 1.0
    if self._inputs["unsupervised_2d"].sum():
      weights["i_con"] = 1.0
    if self._inputs["supervised_2d"].sum():
      weights["dgram_cce"] = 1.0
    self.opt["weights"].update(weights)

    # prep model
    self._prep_model(**kwargs)

  def _prep_hallucination(self,
    length=100,
    copies=1, 
    **kwargs):

    # parsing inputs
    if isinstance(length,int): length = [length]

    # prep model
    self._prep_contigs([str(L) for L in length], copies=copies, **kwargs)

  def _prep_fixbb(self, 
    pdb_filename, 
    chain="A",
    copies=1,
    #rm_template=False,
    #rm_template_seq=True,
    #rm_template_sc=True,
    #rm_template_ic=False,
    fix_pos=None,
    ignore_missing=True,
    parse_as_homooligomer=False, 
    **kwargs):
    
    # parsing inputs
    contigs = chain.split(",") if isinstance(chain,str) else chain

    # prep model
    self._prep_contigs(contigs, pdb_filename,
      fix_pos=fix_pos, ignore_missing=ignore_missing, 
      copies=copies, parse_as_homooligomer=parse_as_homooligomer, **kwargs)

  def _prep_partial(self,
    pdb_filename,
    chain="A",
    length=None,
    copies=1, 
    parse_as_homooligomer=False,
    pos=None, 
    fix_pos=None,
    #use_sidechains=False,
    #atoms_to_exclude=None,
    #rm_template=False,
    #rm_template_seq=False,
    #rm_template_sc=False,
    #rm_template_ic=False,
    ignore_missing=True,
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
      copies=copies, parse_as_homooligomer=parse_as_homooligomer**kwargs)


  def _prep_binder(self,
    pdb_filename,
    target_chain="A",
    binder_len=50,                                         
    binder_chain=None,

    #rm_target = False,
    #rm_target_seq = False,
    #rm_target_sc = False,
    #rm_binder=True,
    #rm_binder_seq=True,
    #rm_binder_sc=True,
    #rm_template_ic=False,
    #hotspot=None,

    ignore_missing=True, **kwargs):
  
    # parsing inputs
    if isinstance(target_chain,str): target_chain = target_chain.split(",")
    if isinstance(binder_chain,str): binder_chain = binder_chain.split(",")
    if binder_chain is None: binder_chain = []

    # define contigs
    fix_pos = []
    contigs = []
    for a in target_chain:
      contigs.append(a)
      fix_pos.append(a)
    for b in binder_chain:
      contigs.append(b)
    if len(binder_chain) == 0:
      contigs.append(str(binder_len))
      redesign = False
    else:
      redesign = True

    self._prep_contigs(contigs, pdb_filename,
      fix_pos=",".join(fix_pos), ignore_missing=ignore_missing, **kwargs)

    self._target_len = sum(self._lengths[:len(target_chain)])
    self._binder_len = sum(self._lengths) - self._target_len

    self._inputs["supervised"][:self._target_len] = False
    self._inputs["supervised_2d"][:self._target_len,:self._target_len] = False

    # adjust weights
    if redesign:
      self.set_weights(i_con=0, con=0, plddt=0, dgram_cce=1, set_defaults=True)
    else:
      self.set_weights(i_con=1, con=0, plddt=0.1, dgram_cce=0, set_defaults=True)

#######################
# utils
#######################

def prep_pdb(pdb_filename, chain=None,
             offsets=None, lengths=None,
             ignore_missing=False, offset_index=False):
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
    pdb_str = pdb_to_string(pdb_filename, chains=chain, models=[1])
    protein_obj = protein.from_pdb_string(pdb_str, chain_id=chain)
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

      batch = {"aatype":scatter(batch["aatype"],-1),
               "all_atom_positions":scatter(batch["all_atom_positions"]),
               "all_atom_mask":scatter(batch["all_atom_mask"]),
               "residue_index":scatter(batch["residue_index"],-1)}
      
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

def prep_input_features(L, N=1, T=1, one_hot_msa=True):
  '''
  given [L]ength, [N]umber of sequences and number of [T]emplates
  return dictionary of blank features
  '''
  inputs = {'aatype': np.zeros(L,int),
            'target_feat': np.zeros((L,20)),
            'deletion_matrix': np.zeros((N,L)),
            'msa_mask': np.ones((N,L)),
            'seq_mask': np.ones(L),                        
            'residue_index': np.arange(L),
            'asym_id': np.zeros(L,int),
            'sym_id': np.zeros(L,int),
            'entity_id': np.zeros(L,int),

            # for template inputs
            'template_aatype': np.zeros((T,L),int),
            'template_all_atom_mask': np.zeros((T,L,37)),
            'template_all_atom_positions': np.zeros((T,L,37,3)),
            'template_mask': np.zeros(T)
            }
  inputs['msa'] = np.zeros((N,L,22)) if one_hot_msa else np.zeros((N,L),int)
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