import jax
import jax.numpy as jnp
import numpy as np
import re

from colabdesign.af.alphafold.data import pipeline, prep_inputs
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

  def _prep_features(self, num_res, num_seq=None, num_templates=1):
    '''process features'''
    if num_seq is None: num_seq = self._num
    return prep_input_features(L=num_res, N=num_seq, T=num_templates)

  def _prep_fixbb(self, pdb_filename, chain=None,
                  copies=1, repeat=False, homooligomer=False,
                  rm_template=False,
                  rm_template_seq=True,
                  rm_template_sc=True,
                  rm_template_ic=False,
                  fix_pos=None, ignore_missing=True, **kwargs):
    '''
    prep inputs for fixed backbone design
    ---------------------------------------------------
    if copies > 1:
      -homooligomer=True - input pdb chains are parsed as homo-oligomeric units
      -repeat=True       - tie the repeating sequence within single chain
    -rm_template_seq     - if template is defined, remove information about template sequence
    -fix_pos="1,2-10"    - specify which positions to keep fixed in the sequence
                           note: supervised loss is applied to all positions, use "partial" 
                           protocol to apply supervised loss to only subset of positions
    -ignore_missing=True - skip positions that have missing density (no CA coordinate)
    ---------------------------------------------------
    '''    
    # prep features
    self._pdb = prep_pdb(pdb_filename, chain=chain, ignore_missing=ignore_missing,
                         offsets=kwargs.pop("pdb_offsets",None),
                         lengths=kwargs.pop("pdb_lengths",None))

    self._len = self._pdb["residue_index"].shape[0]
    self._lengths = [self._len]

    # feat dims
    num_seq = self._num
    res_idx = self._pdb["residue_index"]
    
    # get [pos]itions of interests    
    if fix_pos is not None and fix_pos != "":
      self._pos_info = prep_pos(fix_pos, **self._pdb["idx"])
      self.opt["fix_pos"] = self._pos_info["pos"]

    if homooligomer and chain is not None and copies == 1:
      copies = len(chain.split(","))
      
    # repeat/homo-oligomeric support
    if copies > 1:

      if repeat or homooligomer:
        self._len = self._len // copies
        if "fix_pos" in self.opt:
          self.opt["fix_pos"] = self.opt["fix_pos"][self.opt["fix_pos"] < self._len]

      if repeat:
        self._lengths = [self._len * copies]
        block_diag = False

      else:
        self._lengths = [self._len] * copies
        block_diag = not self._args["use_multimer"]

        res_idx = repeat_idx(res_idx[:self._len], copies)
        num_seq = (self._num * copies + 1) if block_diag else self._num
        self.opt["weights"].update({"i_pae":0.0, "i_con":0.0})

      self._args.update({"copies":copies, "repeat":repeat, "homooligomer":homooligomer, "block_diag":block_diag})
      homooligomer = not repeat
    else:
      self._lengths = self._pdb["lengths"]

    # configure input features
    self._inputs = self._prep_features(num_res=sum(self._lengths), num_seq=num_seq)
    self._inputs["residue_index"] = res_idx
    self._inputs["batch"] = make_fixed_size(self._pdb["batch"], num_res=sum(self._lengths))
    self._inputs.update(get_multi_id(self._lengths, homooligomer=homooligomer))

    # configure options/weights
    self.opt["weights"].update({"dgram_cce":1.0, "rmsd":0.0, "fape":0.0, "con":0.0})
    self._wt_aatype = self._inputs["batch"]["aatype"][:self._len]

    # configure template [opt]ions
    rm,L = {},sum(self._lengths)
    for n,x in {"rm_template":    rm_template,
                "rm_template_seq":rm_template_seq,
                "rm_template_sc": rm_template_sc}.items():
      rm[n] = np.full(L,False)
      if isinstance(x,str):
        rm[n][prep_pos(x,**self._pdb["idx"])["pos"]] = True
      else:
        rm[n][:] = x
    self.opt["template"]["rm_ic"] = rm_template_ic
    self._inputs.update(rm)
  
    self._prep_model(**kwargs)
    
  def _prep_hallucination(self, length=100, copies=1, repeat=False, **kwargs):
    '''
    prep inputs for hallucination
    ---------------------------------------------------
    if copies > 1:
      -repeat=True - tie the repeating sequence within single chain
    ---------------------------------------------------
    '''
    
    # define num copies (for repeats/ homo-oligomers)
    if not repeat and copies > 1 and not self._args["use_multimer"]:
      (num_seq, block_diag) = (self._num * copies + 1, True)
    else:
      (num_seq, block_diag) = (self._num, False)
    
    self._args.update({"repeat":repeat,"block_diag":block_diag,"copies":copies})
      
    # prep features
    self._len = length
    
    # set weights
    self.opt["weights"].update({"con":1.0})
    if copies > 1:
      if repeat:
        offset = 1
        self._lengths = [self._len * copies]
        self._args["repeat"] = True
      else:
        offset = 50
        self._lengths = [self._len] * copies
        self.opt["weights"].update({"i_pae":0.0, "i_con":1.0})
        self._args["homooligomer"] = True
      res_idx = repeat_idx(np.arange(length), copies, offset=offset)
    else:
      self._lengths = [self._len]
      res_idx = np.arange(length)
    
    # configure input features
    self._inputs = self._prep_features(num_res=sum(self._lengths), num_seq=num_seq)
    self._inputs["residue_index"] = res_idx
    self._inputs.update(get_multi_id(self._lengths, homooligomer=True))

    self._prep_model(**kwargs)

  def _prep_binder(self, pdb_filename,
                   target_chain="A", binder_len=50,                                         
                   rm_target = False,
                   rm_target_seq = False,
                   rm_target_sc = False,
                   
                   # if binder_chain is defined
                   binder_chain=None,
                   rm_binder=True,
                   rm_binder_seq=True,
                   rm_binder_sc=True,
                   rm_template_ic=False,
                                      
                   hotspot=None, ignore_missing=True, **kwargs):
    '''
    prep inputs for binder design
    ---------------------------------------------------
    -binder_len = length of binder to hallucinate (option ignored if binder_chain is defined)
    -binder_chain = chain of binder to redesign
    -use_binder_template = use binder coordinates as template input
    -rm_template_ic = use target and binder coordinates as seperate template inputs
    -hotspot = define position/hotspots on target
    -rm_[binder/target]_seq = remove sequence info from template
    -rm_[binder/target]_sc  = remove sidechain info from template
    -ignore_missing=True - skip positions that have missing density (no CA coordinate)
    ---------------------------------------------------
    '''
    redesign = binder_chain is not None
    rm_binder = not kwargs.pop("use_binder_template", not rm_binder)
    
    self._args.update({"redesign":redesign})

    # get pdb info
    target_chain = kwargs.pop("chain",target_chain) # backward comp
    chains = f"{target_chain},{binder_chain}" if redesign else target_chain
    im = [True] * len(target_chain.split(",")) 
    if redesign: im += [ignore_missing] * len(binder_chain.split(","))

    self._pdb = prep_pdb(pdb_filename, chain=chains, ignore_missing=im)
    res_idx = self._pdb["residue_index"]

    if redesign:
      self._target_len = sum([(self._pdb["idx"]["chain"] == c).sum() for c in target_chain.split(",")])
      self._binder_len = sum([(self._pdb["idx"]["chain"] == c).sum() for c in binder_chain.split(",")])
    else:
      self._target_len = self._pdb["residue_index"].shape[0]
      self._binder_len = binder_len
      res_idx = np.append(res_idx, res_idx[-1] + np.arange(binder_len) + 50)
    
    self._len = self._binder_len
    self._lengths = [self._target_len, self._binder_len]

    # gather hotspot info
    if hotspot is not None:
      self.opt["hotspot"] = prep_pos(hotspot, **self._pdb["idx"])["pos"]

    if redesign:
      # binder redesign
      self._wt_aatype = self._pdb["batch"]["aatype"][self._target_len:]
      self.opt["weights"].update({"dgram_cce":1.0, "rmsd":0.0, "fape":0.0,
                                  "con":0.0, "i_con":0.0, "i_pae":0.0})
    else:
      # binder hallucination
      self._pdb["batch"] = make_fixed_size(self._pdb["batch"], num_res=sum(self._lengths))
      self.opt["weights"].update({"plddt":0.1, "con":0.0, "i_con":1.0, "i_pae":0.0})

    # configure input features
    self._inputs = self._prep_features(num_res=sum(self._lengths), num_seq=1)
    self._inputs["residue_index"] = res_idx
    self._inputs["batch"] = self._pdb["batch"]
    self._inputs.update(get_multi_id(self._lengths))

    # configure template rm masks
    (T,L,rm) = (self._lengths[0],sum(self._lengths),{})
    rm_opt = {
              "rm_template":    {"target":rm_target,    "binder":rm_binder},
              "rm_template_seq":{"target":rm_target_seq,"binder":rm_binder_seq},
              "rm_template_sc": {"target":rm_target_sc, "binder":rm_binder_sc}
             }
    for n,x in rm_opt.items():
      rm[n] = np.full(L,False)
      for m,y in x.items():
        if isinstance(y,str):
          rm[n][prep_pos(y,**self._pdb["idx"])["pos"]] = True
        else:
          if m == "target": rm[n][:T] = y
          if m == "binder": rm[n][T:] = y
        
    # set template [opt]ions
    self.opt["template"]["rm_ic"] = rm_template_ic
    self._inputs.update(rm)

    self._prep_model(**kwargs)

  def _prep_partial(self, pdb_filename, chain=None, length=None,
                    copies=1, repeat=False, homooligomer=False,
                    pos=None, fix_pos=None, use_sidechains=False, atoms_to_exclude=None,
                    rm_template=False,
                    rm_template_seq=False,
                    rm_template_sc=False,
                    rm_template_ic=False, 
                    ignore_missing=True, **kwargs):
    '''
    prep input for partial hallucination
    ---------------------------------------------------
    -length=100 - total length of protein (if different from input PDB)
    -pos="1,2-10" - specify which positions to apply supervised loss to
    -use_sidechains=True - add a sidechain supervised loss to the specified positions
      -atoms_to_exclude=["N","C","O"] (for sc_rmsd loss, specify which atoms to exclude)
    -rm_template_seq - if template is defined, remove information about template sequence
    -ignore_missing=True - skip positions that have missing density (no CA coordinate)
    ---------------------------------------------------    
    '''    
    # prep features
    self._pdb = prep_pdb(pdb_filename, chain=chain, ignore_missing=ignore_missing,
                   offsets=kwargs.pop("pdb_offsets",None),
                   lengths=kwargs.pop("pdb_lengths",None))

    self._pdb["len"] = sum(self._pdb["lengths"])

    self._len = self._pdb["len"] if length is None else length
    self._lengths = [self._len]

    # feat dims
    num_seq = self._num
    res_idx = np.arange(self._len)
    
    # get [pos]itions of interests
    if pos is None:
      self.opt["pos"] = self._pdb["pos"] = np.arange(self._pdb["len"])
      self._pos_info = {"length":np.array([self._pdb["len"]]), "pos":self._pdb["pos"]}    
    else:
      self._pos_info = prep_pos(pos, **self._pdb["idx"])
      self.opt["pos"] = self._pdb["pos"] = self._pos_info["pos"]

    if homooligomer and chain is not None and copies == 1:
      copies = len(chain.split(","))

    # repeat/homo-oligomeric support
    if copies > 1:
      
      if repeat or homooligomer:
        self._len = self._len // copies
        self._pdb["len"] = self._pdb["len"] // copies
        self.opt["pos"] = self._pdb["pos"][self._pdb["pos"] < self._pdb["len"]]

        # repeat positions across copies
        self._pdb["pos"] = repeat_pos(self.opt["pos"], copies, self._pdb["len"])

      if repeat:
        self._lengths = [self._len * copies]
        block_diag = False

      else:
        self._lengths = [self._len] * copies
        block_diag = not self._args["use_multimer"]

        num_seq = (self._num * copies + 1) if block_diag else self._num
        res_idx = repeat_idx(np.arange(self._len), copies)

        self.opt["weights"].update({"i_pae":0.0, "i_con":1.0})

      self._args.update({"copies":copies, "repeat":repeat, "homooligomer":homooligomer, "block_diag":block_diag})
      homooligomer = not repeat

    # configure input features
    self._inputs = self._prep_features(num_res=sum(self._lengths), num_seq=num_seq)
    self._inputs["residue_index"] = res_idx
    self._inputs["batch"] = jax.tree_map(lambda x:x[self._pdb["pos"]], self._pdb["batch"])     
    self._inputs.update(get_multi_id(self._lengths, homooligomer=homooligomer))

    # configure options/weights
    self.opt["weights"].update({"dgram_cce":1.0, "rmsd":0.0, "fape":0.0, "con":1.0}) 
    self._wt_aatype = self._pdb["batch"]["aatype"][self.opt["pos"]]

    # configure sidechains
    self._args["use_sidechains"] = use_sidechains
    if use_sidechains:
      self._sc = {"batch":prep_inputs.make_atom14_positions(self._inputs["batch"]),
                  "pos":get_sc_pos(self._wt_aatype, atoms_to_exclude)}
      self.opt["weights"].update({"sc_rmsd":0.1, "sc_fape":0.1})
      self.opt["fix_pos"] = np.arange(self.opt["pos"].shape[0])      
      self._wt_aatype_sub = self._wt_aatype
      
    elif fix_pos is not None and fix_pos != "":
      sub_fix_pos = []
      sub_i = []
      pos = self.opt["pos"].tolist()
      for i in prep_pos(fix_pos, **self._pdb["idx"])["pos"]:
        if i in pos:
          sub_i.append(i)
          sub_fix_pos.append(pos.index(i))
      self.opt["fix_pos"] = np.array(sub_fix_pos)
      self._wt_aatype_sub = self._pdb["batch"]["aatype"][sub_i]
      
    elif kwargs.pop("fix_seq",False):
      self.opt["fix_pos"] = np.arange(self.opt["pos"].shape[0])
      self._wt_aatype_sub = self._wt_aatype

    self.opt["template"].update({"rm_ic":rm_template_ic})
    self._inputs.update({"rm_template":     rm_template,
                         "rm_template_seq": rm_template_seq,
                         "rm_template_sc":  rm_template_sc})
  
    self._prep_model(**kwargs)

#######################
# utils
#######################
def repeat_idx(idx, copies=1, offset=50):
  idx_offset = np.repeat(np.cumsum([0]+[idx[-1]+offset]*(copies-1)),len(idx))
  return np.tile(idx,copies) + idx_offset

def repeat_pos(pos, copies, length):
  return (np.repeat(pos,copies).reshape(-1,copies) + np.arange(copies) * length).T.flatten()

def prep_pdb(pdb_filename, chain=None,
             offsets=None, lengths=None,
             ignore_missing=False):
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

  if isinstance(chain,str) and "," in chain:
    chains = chain.split(",")
  elif not isinstance(chain,list):
    chains = [chain]

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
      r = batch["all_atom_mask"][:,0] == 1
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
    
    last = residue_index[-1] + 50
    o.append({"batch":batch,
              "residue_index": residue_index,
              "cb_feat":cb_feat})
    
    residue_idx.append(batch.pop("residue_index"))
    chain_idx.append([chain] * len(residue_idx[-1]))
    full_lengths.append(len(residue_index))

  # concatenate chains
  o = jax.tree_util.tree_map(lambda *x:np.concatenate(x,0),*o)
  
  # save original residue and chain index
  o["idx"] = {"residue":np.concatenate(residue_idx), "chain":np.concatenate(chain_idx)}
  o["lengths"] = full_lengths
  return o

def make_fixed_size(feat, num_res, num_seq=1, num_templates=1):
  '''pad input features'''
  shape_schema = {k:v for k,v in config.CONFIG.data.eval.feat.items()}

  pad_size_map = {
      shape_placeholders.NUM_RES: num_res,
      shape_placeholders.NUM_MSA_SEQ: num_seq,
      shape_placeholders.NUM_EXTRA_SEQ: 1,
      shape_placeholders.NUM_TEMPLATES: num_templates
  }  
  for k,v in feat.items():
    if k == "batch":
      feat[k] = make_fixed_size(v, num_res)
    else:
      shape = list(v.shape)
      schema = shape_schema[k]
      assert len(shape) == len(schema), (
          f'Rank mismatch between shape and shape schema for {k}: '
          f'{shape} vs {schema}')
      pad_size = [pad_size_map.get(s2, None) or s1 for (s1, s2) in zip(shape, schema)]
      padding = [(0, p - v.shape[i]) for i, p in enumerate(pad_size)]
      feat[k] = np.pad(v, padding)
  return feat

def get_sc_pos(aa_ident, atoms_to_exclude=None):
  '''get sidechain indices/weights for all_atom14_positions'''

  # decide what atoms to exclude for each residue type
  a2e = {}
  for r in resname_to_idx:
    if isinstance(atoms_to_exclude,dict):
      a2e[r] = atoms_to_exclude.get(r,atoms_to_exclude.get("ALL",["N","C","O"]))
    else:
      a2e[r] = ["N","C","O"] if atoms_to_exclude is None else atoms_to_exclude

  # collect atom indices
  pos,pos_alt = [],[]
  N,N_non_amb = [],[]
  for n,a in enumerate(aa_ident):
    aa = idx_to_resname[a]
    atoms = set(residue_constants.residue_atoms[aa])
    atoms14 = residue_constants.restype_name_to_atom14_names[aa]
    swaps = residue_constants.residue_atom_renaming_swaps.get(aa,{})
    swaps.update({v:k for k,v in swaps.items()})
    for atom in atoms.difference(a2e[aa]):
      pos.append(n * 14 + atoms14.index(atom))
      if atom in swaps:
        pos_alt.append(n * 14 + atoms14.index(swaps[atom]))
      else:
        pos_alt.append(pos[-1])
        N_non_amb.append(n)
      N.append(n)

  pos, pos_alt = np.asarray(pos), np.asarray(pos_alt)
  non_amb = pos == pos_alt
  N, N_non_amb = np.asarray(N), np.asarray(N_non_amb)
  w = np.array([1/(n == N).sum() for n in N])
  w_na = np.array([1/(n == N_non_amb).sum() for n in N_non_amb])
  w, w_na = w/w.sum(), w_na/w_na.sum()
  return {"pos":pos, "pos_alt":pos_alt, "non_amb":non_amb,
          "weight":w, "weight_non_amb":w_na[:,None]}

def prep_input_features(L, N=1, T=1, eN=1):
  '''
  given [L]ength, [N]umber of sequences and number of [T]emplates
  return dictionary of blank features
  '''
  inputs = {'aatype': np.zeros(L,int),
            'target_feat': np.zeros((L,20)),
            'msa_feat': np.zeros((N,L,49)),
            # 23 = one_hot -> (20, UNK, GAP, MASK)
            # 1  = has deletion
            # 1  = deletion_value
            # 23 = profile
            # 1  = deletion_mean_value
  
            'seq_mask': np.ones(L),
            'msa_mask': np.ones((N,L)),
            'msa_row_mask': np.ones(N),
            'atom14_atom_exists': np.zeros((L,14)),
            'atom37_atom_exists': np.zeros((L,37)),
            'residx_atom14_to_atom37': np.zeros((L,14),int),
            'residx_atom37_to_atom14': np.zeros((L,37),int),            
            'residue_index': np.arange(L),
            'extra_deletion_value': np.zeros((eN,L)),
            'extra_has_deletion': np.zeros((eN,L)),
            'extra_msa': np.zeros((eN,L),int),
            'extra_msa_mask': np.zeros((eN,L)),
            'extra_msa_row_mask': np.zeros(eN),

            # for template inputs
            'template_aatype': np.zeros((T,L),int),
            'template_all_atom_mask': np.zeros((T,L,37)),
            'template_all_atom_positions': np.zeros((T,L,37,3)),
            'template_mask': np.zeros(T),
            'template_pseudo_beta': np.zeros((T,L,3)),
            'template_pseudo_beta_mask': np.zeros((T,L)),

            # for alphafold-multimer
            'asym_id': np.zeros(L),
            'sym_id': np.zeros(L),
            'entity_id': np.zeros(L),
            'all_atom_positions': np.zeros((N,37,3))}
  return inputs

def get_multi_id(lengths, homooligomer=False):
  '''set info for alphafold-multimer'''
  i = np.concatenate([[n]*l for n,l in enumerate(lengths)])
  if homooligomer:
    return {"asym_id":i, "sym_id":i, "entity_id":np.zeros_like(i)}
  else:
    return {"asym_id":i, "sym_id":i, "entity_id":i}