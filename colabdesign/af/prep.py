import tensorflow as tf
tf.config.set_visible_devices([], 'GPU')

import jax
import jax.numpy as jnp
import numpy as np
import re

from colabdesign.af.alphafold.data import pipeline, prep_inputs
from colabdesign.af.alphafold.common import protein, residue_constants
from colabdesign.af.alphafold.model import all_atom
from colabdesign.af.alphafold.model.tf import shape_placeholders

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

  def _prep_features(self, length, num_templates=1, template_features=None):
    '''process features'''
    return prep_input_features(L=length, N=self._num, T=num_templates,
                               use_templates=self._args["use_templates"])
  
  # prep functions specific to protocol
  def _prep_binder(self, pdb_filename, chain="A",
                   binder_len=50, binder_chain=None,
                   use_binder_template=False, split_templates=False,
                   pos=None, rm_template_seq=True, **kwargs):
    '''
    prep inputs for binder design
    ---------------------------------------------------
    -binder_len = length of binder to hallucinate (option ignored if binder_chain is defined)
    -binder_chain = chain of binder to redesign
    -use_binder_template = use binder coordinates as template input
    -split_templates = use target and binder coordinates as seperate template inputs
    -pos = define position/hotspots on target
    -rm_template_seq = for binder redesign protocol, remove sequence info from binder template
    '''
    
    redesign = binder_chain is not None

    self._args.update({"rm_template_seq":rm_template_seq,
                       "redesign":redesign})

    self.opt["template"]["dropout"] = 0.0 if use_binder_template else 1.0
    num_templates = 1

    # get pdb info
    chains = f"{chain},{binder_chain}" if redesign else chain
    pdb = prep_pdb(pdb_filename, chain=chains)
    if redesign:
      target_len = sum([(pdb["idx"]["chain"] == c).sum() for c in chain.split(",")])
      binder_len = sum([(pdb["idx"]["chain"] == c).sum() for c in binder_chain.split(",")])
      if split_templates: num_templates = 2      
      # get input features
      self._inputs = self._prep_features(target_len + binder_len, num_templates=num_templates)

    else:
      target_len = pdb["residue_index"].shape[0]
      self._inputs = self._prep_features(target_len)
      
    self._inputs["residue_index"][...,:] = pdb["residue_index"]

    # gather hotspot info
    hotspot = kwargs.pop("hotspot", pos)
    if hotspot is not None:
      self.opt["pos"] = prep_pos(hotspot, **pdb["idx"])["pos"]

    if redesign:      
      self._batch = pdb["batch"]
      self._wt_aatype = self._batch["aatype"][target_len:]
      self.opt["weights"].update({"dgram_cce":1.0, "fape":0.0, "rmsd":0.0,
                                  "con":0.0, "i_pae":0.01, "i_con":0.0})      
    else: # binder hallucination            
      # pad inputs
      total_len = target_len + binder_len
      self._inputs = make_fixed_size(self._inputs, self._runner, total_len)
      self._batch = make_fixed_size(pdb["batch"], self._runner, total_len, batch_axis=False)

      # offset residue index for binder
      self._inputs["residue_index"] = self._inputs["residue_index"].copy()
      self._inputs["residue_index"][:,target_len:] = pdb["residue_index"][-1] + np.arange(binder_len) + 50
      for k in ["seq_mask","msa_mask"]: self._inputs[k] = np.ones_like(self._inputs[k])
      self.opt["weights"].update({"con":0.5, "i_pae":0.01, "i_con":0.5})

    self._target_len = target_len
    self._binder_len = self._len = binder_len

    self._opt = copy_dict(self.opt)
    self.restart(**kwargs)

  def _prep_fixbb(self, pdb_filename, chain=None, copies=1, homooligomer=False, 
                  repeat=False, block_diag=False, rm_template_seq=True,
                  pos=None, fix_seq=False, **kwargs):
    '''prep inputs for fixed backbone design'''

    self._args["rm_template_seq"] = rm_template_seq

    # block_diag the msa features
    if block_diag and not repeat and copies > 1:
      self._runner.config.data.eval.max_msa_clusters = self._num * (1 + copies)
    else:
      block_diag = False

    pdb = prep_pdb(pdb_filename, chain=chain)
    self._batch = pdb["batch"]
    self._len = pdb["residue_index"].shape[0]
    self._inputs = self._prep_features(self._len)
    self._copies = copies
    self._args.update({"repeat":repeat, "block_diag":block_diag, "homooligomer":homooligomer})

    # set weights
    self.opt["weights"].update({"dgram_cce":1.0, "rmsd":0.0, "con":0.0, "fape":0.0})

    # update residue index from pdb
    if copies > 1:
      if repeat:
        self._len = self._len // copies
        block_diag = False
      else:
        if homooligomer:
          self._len = self._len // copies
          self._inputs["residue_index"] = pdb["residue_index"][None]
        else:
          self._inputs = make_fixed_size(self._inputs, self._runner, self._len * copies)
          self._batch = make_fixed_size(self._batch, self._runner, self._len * copies, batch_axis=False)
          self._inputs["residue_index"] = repeat_idx(pdb["residue_index"], copies)[None]
          for k in ["seq_mask","msa_mask"]: self._inputs[k] = np.ones_like(self._inputs[k])
    else:
      self._inputs["residue_index"] = pdb["residue_index"][None]

    # fix certain positions
    self.opt["fix_seq"] = fix_seq
    if pos is not None:
      self._pos_info = prep_pos(pos, **pdb["idx"])
      self.opt["pos"] = self._pos_info["pos"]

    self._wt_aatype = self._batch["aatype"][:self._len]
    self._opt = copy_dict(self.opt)
    self.restart(**kwargs)
    
  def _prep_hallucination(self, length=100, copies=1,
                          repeat=False, block_diag=False, **kwargs):
    '''prep inputs for hallucination'''
    
    # block_diag the msa features
    if block_diag and not repeat and copies > 1:
      self._runner.config.data.eval.max_msa_clusters = self._num * (1 + copies)
    else:
      block_diag = False
      
    self._len = length
    self._copies = copies
    self._inputs = self._prep_features(length * copies)
    self._args.update({"block_diag":block_diag, "repeat":repeat})
    
    # set weights
    self.opt["weights"].update({"con":1.0})
    if copies > 1:
      if repeat:
        offset = 1
      else:
        offset = 50
        self.opt["weights"].update({"i_pae":0.01, "i_con":0.1})
      self._inputs["residue_index"] = repeat_idx(np.arange(length), copies, offset=offset)[None]

    self._opt = copy_dict(self.opt)
    self.restart(**kwargs)

  def _prep_partial(self, pdb_filename, chain=None, length=None,
                    pos=None, fix_seq=True, use_sidechains=False, atoms_to_exclude=None,
                    rm_template_seq=False, **kwargs):
    '''prep input for partial hallucination'''    

    self._args["rm_template_seq"] = rm_template_seq

    # prep features
    pdb = prep_pdb(pdb_filename, chain=chain)
    self._len = pdb["residue_index"].shape[0] if length is None else length
    self._inputs = self._prep_features(self._len)    

    # configure options/weights
    self.opt["weights"].update({"dgram_cce":1.0,"con":1.0, "fape":0.0, "rmsd":0.0})
    self.opt["fix_seq"] = fix_seq

    # get [pos]itions of interests
    if pos is None:
      self.opt["pos"] = np.arange(pdb["residue_index"].shape[0])
    else:
      self._pos_info = prep_pos(pos, **pdb["idx"])
      self.opt["pos"] = p = self._pos_info["pos"]
      self._batch = jax.tree_map(lambda x:x[p], pdb["batch"])     
    self._wt_aatype = self._batch["aatype"]

    # configure sidechains
    self._args["use_sidechains"] = kwargs.pop("sidechain", use_sidechains)
    if self._args["use_sidechains"]:
      self._batch.update(prep_inputs.make_atom14_positions(self._batch))
      self._batch["sc_pos"] = get_sc_pos(self._wt_aatype, atoms_to_exclude)
      self.opt["weights"].update({"sc_rmsd":0.1, "sc_fape":0.1})
      self.opt["fix_seq"] = True
  
    self._opt = copy_dict(self.opt)
    self.restart(**kwargs)

#######################
# utils
#######################
def repeat_idx(idx, copies=1, offset=50):
  idx_offset = np.repeat(np.cumsum([0]+[idx[-1]+offset]*(copies-1)),len(idx))
  return np.tile(idx,copies) + idx_offset

def prep_pdb(pdb_filename, chain=None, for_alphafold=True):
  '''extract features from pdb'''

  def add_cb(batch):
    '''add missing CB atoms based on N,CA,C'''
    p,m = batch["all_atom_positions"],batch["all_atom_mask"]
    atom_idx = residue_constants.atom_order
    atoms = {k:p[...,atom_idx[k],:] for k in ["N","CA","C"]}
    cb = atom_idx["CB"]
    cb_atoms = _np_get_cb(**atoms, use_jax=False)
    cb_mask = np.prod([m[...,atom_idx[k]] for k in ["N","CA","C"]],0)
    batch["all_atom_positions"][...,cb,:] = np.where(m[:,cb,None], p[:,cb,:], cb_atoms)
    batch["all_atom_mask"][...,cb] = (m[:,cb] + cb_mask) > 0

  # go through each defined chain
  chains = [None] if chain is None else chain.split(",")
  o,last = [],0
  residue_idx, chain_idx = [],[]
  for chain in chains:
    protein_obj = protein.from_pdb_string(pdb_to_string(pdb_filename), chain_id=chain)
    batch = {'aatype': protein_obj.aatype,
             'all_atom_positions': protein_obj.atom_positions,
             'all_atom_mask': protein_obj.atom_mask}

    add_cb(batch) # add in missing cb (in the case of glycine)

    has_ca = batch["all_atom_mask"][:,0] == 1
    batch = jax.tree_map(lambda x:x[has_ca], batch)
    seq = "".join([order_aa[a] for a in batch["aatype"]])
    residue_index = protein_obj.residue_index[has_ca] + last      
    last = residue_index[-1] + 50
    
    if for_alphafold:
      batch.update(all_atom.atom37_to_frames(**batch))
      template_aatype = residue_constants.sequence_to_onehot(seq, residue_constants.HHBLITS_AA_TO_ID)
      template_features = {"template_aatype":template_aatype,
                           "template_all_atom_masks":batch["all_atom_mask"],
                           "template_all_atom_positions":batch["all_atom_positions"]}
      o.append({"batch":batch,
                "template_features":template_features,
                "residue_index": residue_index})
    else:        
      o.append({"batch":batch,
                "residue_index": residue_index})
    
    residue_idx.append(protein_obj.residue_index[has_ca])
    chain_idx.append([chain] * len(residue_idx[-1]))

  # concatenate chains
  o = jax.tree_util.tree_map(lambda *x:np.concatenate(x,0),*o)
  
  if for_alphafold:
    o["template_features"] = jax.tree_map(lambda x:x[None],o["template_features"])
    o["template_features"]["template_domain_names"] = np.asarray(["None"])

  # save original residue and chain index
  o["idx"] = {"residue":np.concatenate(residue_idx), "chain":np.concatenate(chain_idx)}
  return o

def make_fixed_size(feat, model_runner, length, batch_axis=True):
  '''pad input features'''
  cfg = model_runner.config
  if batch_axis:
    shape_schema = {k:[None]+v for k,v in dict(cfg.data.eval.feat).items()}
  else:
    shape_schema = {k:v for k,v in dict(cfg.data.eval.feat).items()}
  num_msa_seq = cfg.data.eval.max_msa_clusters - cfg.data.eval.max_templates
  pad_size_map = {
      shape_placeholders.NUM_RES: length,
      shape_placeholders.NUM_MSA_SEQ: num_msa_seq,
      shape_placeholders.NUM_EXTRA_SEQ: cfg.data.common.max_extra_msa,
      shape_placeholders.NUM_TEMPLATES: cfg.data.eval.max_templates
  }
  for k, v in feat.items():
    # Don't transfer this to the accelerator.
    if k == 'extra_cluster_assignment':
      continue
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

def prep_input_features(L, N=1, T=1, use_templates=False, eN=1):
  '''
  given [L]ength, [N]umber of sequences and number of [T]emplates
  return dictionary of blank features
  '''
  inputs = {'aatype': np.zeros(L,int),
            'target_feat': np.zeros((L,22)),
            'msa_feat': np.zeros((N,L,49)),
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
            'extra_msa_row_mask': np.zeros(eN)}
                   
  if use_templates:
    inputs.update({'template_aatype': np.zeros((T,L),int),
                   'template_all_atom_masks': np.zeros((T,L,37)),
                   'template_all_atom_positions': np.zeros((T,L,37,3)),
                   'template_mask': np.ones(T),
                   'template_pseudo_beta': np.zeros((T,L,3)),
                   'template_pseudo_beta_mask': np.zeros((T,L))})

  return jax.tree_map(lambda x:x[None], inputs)