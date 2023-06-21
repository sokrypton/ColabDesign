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

  def _prep_fixbb(self, 
                  pdb_filename, 
                  chain="A",
                  copies=1,
                  rm_template=False,
                  rm_template_seq=True,
                  rm_template_sc=True,
                  rm_template_ic=False,
                  fix_pos=None,
                  ignore_missing=True,
                  parse_as_homooligomer=False, 
                  **kwargs):
    '''
    prep inputs for fixed backbone design
    ---------------------------------------------------
    -copies                - number of copies in homooligomeric unit
    -parse_as_homooligomer - input pdb chains are parsed as homo-oligomeric units
    -rm_template_seq       - if template is defined, remove information about template sequence
    -fix_pos="1,2-10"      - specify which positions to keep fixed in the sequence
    -ignore_missing=True   - skip positions that have missing density (no CA coordinate)
    ---------------------------------------------------
    '''    
    
    parse_as_homooligomer = kwargs.pop("homooligomer",parse_as_homooligomer)    
    if isinstance(chain,str): chain = chain.split(",")
    
    # prep features
    self._pdb = prep_pdb(pdb_filename, chain=chain, ignore_missing=ignore_missing)

    if parse_as_homooligomer:
      if copies == 1: copies = len(chain)
      assert len(self._pdb["lengths"]) % copies == 0
      units = len(self._pdb["lengths"]) // copies
      length = self._pdb["lengths"][:units]
      residue_index = self._pdb["residue_index"][:sum(length)]
    else: 
      length = self._pdb["lengths"]
      residue_index = self._pdb["residue_index"]

    # encode chains
    if isinstance(length, int): length = [length]
    self._len = sum(length)
    self._lengths, chain_enc = get_chain_enc(copies, length, residue_index)
    
    # get [pos]itions of interests 
    if fix_pos is not None:
      fix_pos = prep_pos(fix_pos, **self._pdb["idx"])["pos"]
      if parse_as_homooligomer:
        fix_pos = fix_pos[fix_pos < self._len]
      self.opt["fix_pos"] = fix_pos
    
    self._args.update({
      "copies":copies,
      "block_diag":copies > 1 and not self._args["use_multimer"]:
    })

    # configure input features
    self._inputs = self._prep_features(num_res=sum(self._lengths))
    self._inputs["batch"] = make_fixed_size(self._pdb["batch"],
      num_res=sum(self._lengths), num_templates=self._args["num_templates"])
    self._inputs.update(chain_enc)

    # configure options/weights
    self.opt["weights"].update({"dgram_cce":1.0, "rmsd":0.0, "fape":0.0, "con":0.0})
    if len(self._lengths) > 1:
      self.opt["weights"].update({"i_pae":0.0, "i_con":0.0})

    # add sidechain loss
    if kwargs.pop("use_sidechains", self._args["use_sidechains"]):
      self._args["use_sidechains"] = True
      self.opt["weights"].update({"sc_fape":0.1, "sc_chi":0.1, "sc_chi_norm":0.0})

    self._wt_aatype = self._inputs["batch"]["aatype"][:self._len]

    # configure template masks
    if self._args["use_templates"]:
      rm_dict = {}
      rm_opt = {"rm_template":    rm_template,
                "rm_template_seq":rm_template_seq,
                "rm_template_sc": rm_template_sc}
      for n,x in rm_opt.items():
        rm_dict[n] = np.full(sum(self._lengths),False)
        if isinstance(x,str):
          rm_dict[n][prep_pos(x,**self._pdb["idx"])["pos"]] = True
        else:
          rm_dict[n][:] = x
      self._inputs.update(rm_dict)
      self.opt["template"]["rm_ic"] = rm_template_ic
  
    self._prep_model(**kwargs)
    
  def _prep_hallucination(self, length=100, copies=1, **kwargs):
    '''
    prep inputs for hallucination
    '''    
    # define num copies (for homo-oligomers)
    self._args.update({
      "block_diag":copies > 1 and not self._args["use_multimer"],
      "copies":copies})   
    
    # encode chains
    if isinstance(length, int): length = [length]
    self._len = sum(length)
    self._lengths, chain_enc = get_chain_enc(copies, length)
    
    # set weights
    self.opt["weights"].update({"con":1.0})
    if copies > 1:
      self.opt["weights"].update({"i_pae":0.0, "i_con":1.0})
    
    # configure input features
    self._inputs = self._prep_features(num_res=sum(self._lengths))
    self._inputs.update(chain_enc)

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
    
    # get pdb info
    target_chain = kwargs.pop("chain",target_chain) # backward comp
    if isinstance(target_chain,str): target_chain = target_chain.split(",")
    if isinstance(binder_chain,str): binder_chain = binder_chain.split(",")
    
    # decide how to parse chains
    if binder_chain is None or len(binder_chain) == 0:
      # binder hallucination
      redesign = False
      chains = target_chain
      ignore_missing = [True] * len(target_chain)
    else:
      # binder redesign
      redesign = True
      chains = target_chain + binder_chain
      ignore_missing = [True] * len(target_chain) + [ignore_missing] * len(binder_chain)
    self._args.update({"redesign":redesign})

    # parse pdb
    self._pdb = prep_pdb(pdb_filename, chain=chains, ignore_missing=ignore_missing)
    target_len = [(self._pdb["idx"]["chain"] == c).sum() for c in target_chain]

    # get lengths
    if redesign: 
      binder_len = [(self._pdb["idx"]["chain"] == c).sum() for c in binder_chain]
    elif not isinstance(binder_len,list): 
      binder_len = [binder_len]
    
    self._target_len = sum(target_len)
    self._binder_len = self._len = sum(binder_len)
    self._lengths = target_len + binder_len

    # chain encoding
    i = np.concatenate([[n]*l for n,l in enumerate(self._lengths)])
    chain_enc = {"asym_id":i, "sym_id":i, "entity_id":i}

    # gather hotspot info
    if hotspot is not None:
      self.opt["hotspot"] = prep_pos(hotspot, **self._pdb["idx"])["pos"]

    if redesign:
      self._wt_aatype = self._pdb["batch"]["aatype"][self._target_len:]
      self.opt["weights"].update({"dgram_cce":1.0, "rmsd":0.0, "fape":0.0,
                                  "con":0.0, "i_con":0.0, "i_pae":0.0})
      # add sidechain loss
      if kwargs.pop("use_sidechains", self._args["use_sidechains"]):
        self._args["use_sidechains"] = True
        self.opt["weights"].update({"sc_fape":0.1, "sc_chi":0.1, "sc_chi_norm":0.0})

    else:
      self._pdb["batch"] = make_fixed_size(self._pdb["batch"],
        num_res=sum(self._lengths), num_templates=self._args["num_templates"])
      self.opt["weights"].update({"plddt":0.1, "con":0.0, "i_con":1.0, "i_pae":0.0})

    # configure input features
    self._inputs = self._prep_features(num_res=sum(self._lengths))
    self._inputs["batch"] = self._pdb["batch"]
    self._inputs.update(chain_enc)

    # configure residue index
    res_idx = self._pdb["residue_index"]
    if not redesign:
      binder_res_idx = np.concatenate([np.arange(L) for L in binder_len])
      res_idx = np.append(res_idx, binder_res_idx)
    self._inputs["residue_index"] = res_idx

    # configure template masks
    rm_binder = not kwargs.pop("use_binder_template", not rm_binder) # backward comp
    rm_dict = {}
    rm_opt = {
              "rm_template":    {"target":rm_target,    "binder":rm_binder},
              "rm_template_seq":{"target":rm_target_seq,"binder":rm_binder_seq},
              "rm_template_sc": {"target":rm_target_sc, "binder":rm_binder_sc}
             }
    for n,x in rm_opt.items():
      rm_dict[n] = np.full(sum(self._lengths), False)
      for m,y in x.items():
        if isinstance(y,str):
          rm_dict[n][prep_pos(y,**self._pdb["idx"])["pos"]] = True
        else:
          if m == "target": rm_dict[n][:self._target_len] = y
          if m == "binder": rm_dict[n][self._target_len:] = y
    self._inputs.update(rm_dict)        
    self.opt["template"]["rm_ic"] = rm_template_ic

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
    
    residue_idx.append(batch.pop("residue_index"))
    chain_idx.append([chain] * len(residue_idx[-1]))
    full_lengths.append(len(residue_index))

  # concatenate chains
  o = jax.tree_util.tree_map(lambda *x:np.concatenate(x,0),*o)
  
  # save original residue and chain index
  o["idx"] = {"residue":np.concatenate(residue_idx), "chain":np.concatenate(chain_idx)}
  o["lengths"] = full_lengths
  return o

def smarter_pad(arr, pad_values, pad_value=0):
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

def make_fixed_size(feat, num_res, num_seq=1, num_templates=1, skip_batch=False):
  '''pad input features'''
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
      feat[k] = smarter_pad(v, padding)
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
    for entity_id, (L,R) in enumerate(zip(lengths,residue_idxs)):
      chain_enc["asym_id"].append(np.full(L,asym_id))
      chain_enc["sym_id"].append(np.full(L,sym_id))
      chain_enc["entity_id"].append(np.full(L,entity_id))
      chain_lengths.append(residue_idx)
      asym_id += 1
    chain_enc["residue_index"].append(R)
  for k in chain_enc.keys():
    chain_enc[k] = np.concatenate(chain_enc[k])
  return chain_lengths, chain_enc