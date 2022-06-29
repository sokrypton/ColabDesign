import tensorflow as tf
tf.config.set_visible_devices([], 'GPU')
import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.af.alphafold.data import pipeline, prep_inputs
from colabdesign.af.alphafold.common import protein, residue_constants
from colabdesign.af.alphafold.model import all_atom
from colabdesign.af.alphafold.model.tf import shape_placeholders

from colabdesign.af.misc import _np_get_cb

ORDER_RESTYPE = {v: k for k, v in residue_constants.restype_order.items()}

#################################################
# AF_PREP - input prep functions
#################################################
class _af_prep:
  def repeat_idx(self, idx, copies=1, offset=50):
    idx_offset = np.repeat(np.cumsum([0]+[idx[-1]+offset]*(copies-1)),len(idx))
    return np.tile(idx,copies) + idx_offset

  def _prep_features(self, length, num_templates=1, template_features=None):
    '''process features'''
    sequence = "A" * length
    feature_dict = {
        **pipeline.make_sequence_features(sequence=sequence, description="none", num_res=length),
        **pipeline.make_msa_features(msas=[length*[sequence]], deletion_matrices=[self._num*[[0]*length]])
    }
    if self._args["use_templates"]:
      # reconfigure model
      self._runner.config.data.eval.max_templates = num_templates
      self._runner.config.data.eval.max_msa_clusters = self._num + num_templates
      if template_features is None:
        feature_dict.update({'template_aatype': np.zeros([num_templates,length,22]),
                             'template_all_atom_masks': np.zeros([num_templates,length,37]),
                             'template_all_atom_positions': np.zeros([num_templates,length,37,3]),
                             'template_domain_names': np.zeros(num_templates)})
      else:
        feature_dict.update(template_features)
      
    inputs = self._runner.process_features(feature_dict, random_seed=0)
    if self._num > 1:
      inputs["msa_row_mask"] = jnp.ones_like(inputs["msa_row_mask"])
      inputs["msa_mask"] = jnp.ones_like(inputs["msa_mask"])
    return inputs
  
  # prep functions specific to protocol
  def _prep_binder(self, pdb_filename, chain="A",
                   binder_len=50, binder_chain=None,
                   use_binder_template=False, split_templates=False,
                   hotspot=None, **kwargs):
    '''
    ---------------------------------------------------
    prep inputs for binder design
    ---------------------------------------------------
    -binder_len = length of binder to hallucinate (option ignored if binder_chain is defined)
    -binder_chain = chain of binder to redesign
    -use_binder_template = use binder coordinates as template input
    -split_templates = use target and binder coordinates as seperate template inputs
    -hotspot = define position on target
    '''
    
    self._redesign = binder_chain is not None
    self._copies = 1
    self._opt["template"]["dropout"] = 0.0 if use_binder_template else 1.0
    num_templates = 1

    # get pdb info
    chains = f"{chain},{binder_chain}" if self._redesign else chain
    pdb = prep_pdb(pdb_filename, chain=chains)
    if self._redesign:
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
    if hotspot is None:
      self._hotspot = None
    else:
      self._hotspot = prep_pos(hotspot, **pdb["idx"])

    if self._redesign:      
      self._batch = pdb["batch"]
      self._wt_aatype = self._batch["aatype"][target_len:]
      self._opt["weights"].update({"dgram_cce":1.0, "fape":0.0, "rmsd":0.0,
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
      self._opt["weights"].update({"con":0.5, "i_pae":0.01, "i_con":0.5})

    self._target_len = target_len
    self._binder_len = self._len = binder_len

    self.restart(set_defaults=True, **kwargs)

  def _prep_fixbb(self, pdb_filename, chain=None, copies=1, 
                  homooligomer=False, repeat=False, block_diag=False, **kwargs):
    '''prep inputs for fixed backbone design'''

    # block_diag the msa features
    if block_diag and not repeat and copies > 1:
      self._runner.config.data.eval.max_msa_clusters = self._num * (1 + copies)
    else:
      block_diag = False

    pdb = prep_pdb(pdb_filename, chain=chain)
    self._batch = pdb["batch"]
    self._wt_aatype = self._batch["aatype"]
    self._len = pdb["residue_index"].shape[0]
    self._inputs = self._prep_features(self._len)
    self._copies = copies
    self._args.update({"repeat":repeat,"block_diag":block_diag})
    
    # set weights
    self._opt["weights"].update({"dgram_cce":1.0, "rmsd":0.0, "con":0.0, "fape":0.0})

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
          self._inputs["residue_index"] = self.repeat_idx(pdb["residue_index"], copies)[None]
          for k in ["seq_mask","msa_mask"]: self._inputs[k] = np.ones_like(self._inputs[k])
        self._opt["weights"].update({"i_pae":0.01, "i_con":0.0})
    else:
      self._inputs["residue_index"] = pdb["residue_index"][None]

    self.restart(set_defaults=True, **kwargs)
    
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
    self._opt["weights"].update({"con":1.0})
    if copies > 1:
      if repeat:
        offset = 1
      else:
        offset = 50
        self._opt["weights"].update({"i_pae":0.01, "i_con":0.1})
      self._inputs["residue_index"] = self.repeat_idx(np.arange(length), copies, offset=offset)[None]

    self.restart(set_defaults=True, **kwargs)

  def _prep_partial(self, pdb_filename, chain=None, pos=None, length=None,
                    fix_seq=False, use_sidechains=False, use_6D=False, **kwargs):
    '''prep input for partial hallucination'''
    
    if "sidechain" in kwargs:
      use_sidechains = kwargs.pop("sidechain")
    self._args["use_sidechains"] = use_sidechains
    if use_sidechains: fix_seq = True

    self._copies = 1
    
    # get [pos]itions of interests
    pdb = prep_pdb(pdb_filename, chain=chain)
    pos = np.arange(pdb["residue_index"].shape[0]) if pos is None else prep_pos(pos, **pdb["idx"])
        
    self._opt.update({"pos":pos, "fix_seq":fix_seq})
    self._batch = jax.tree_map(lambda x:x[pos], pdb["batch"])
    self._wt_aatype = self._batch["aatype"]

    if use_sidechains:
      self._batch.update(prep_inputs.make_atom14_positions(self._batch))

    self._len = pdb["residue_index"].shape[0] if length is None else length
    self._inputs = self._prep_features(self._len)
    
    weights = {"dgram_cce":1.0,"con":1.0, "fape":0.0, "rmsd":0.0}
    if use_sidechains:
      weights.update({"sc_rmsd":0.0, "sc_fape":0.0})
      
    self._opt["weights"].update(weights)
    self.restart(set_defaults=True, **kwargs)

#######################
# utils
#######################
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
    seq = "".join([ORDER_RESTYPE[a] for a in batch["aatype"]])
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

def prep_pos(pos, residue, chain=None):
  if chain is None: chain = np.full(len(residue),None)
  residue_set,chain_set = [],[]
  for idx in pos.split(","):
    i,j = idx.split("-") if "-" in idx else (idx, None)

    # if chain defined
    if i[0].isalpha():
      c,i = i[0], int(i[1:])
    else:
      c,i = chain[0], int(i)

    if j is None:
      residue_set.append(i)
      chain_set.append(c)
    else:
      j = int(j[1:] if j[0].isalpha() else j)
      residue_set += list(range(i,j+1))
      chain_set += [c] * (j-i+1)

  # get positions
  pos_array = []
  for i,c in zip(residue_set, chain_set):
    idx = np.where((residue == i) & (chain == c))[0]
    if idx.shape[0] > 0: pos_array.append(idx[0])

  return np.asarray(pos_array)

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
    padding = [(0, p - tf.shape(v)[i]) for i, p in enumerate(pad_size)]
    if padding:
      feat[k] = tf.pad(v, padding, name=f'pad_to_fixed_{k}')
      feat[k].set_shape(pad_size)
  return {k:np.asarray(v) for k,v in feat.items()}

MODRES = {'MSE':'MET','MLY':'LYS','FME':'MET','HYP':'PRO',
          'TPO':'THR','CSO':'CYS','SEP':'SER','M3L':'LYS',
          'HSK':'HIS','SAC':'SER','PCA':'GLU','DAL':'ALA',
          'CME':'CYS','CSD':'CYS','OCS':'CYS','DPR':'PRO',
          'B3K':'LYS','ALY':'LYS','YCM':'CYS','MLZ':'LYS',
          '4BF':'TYR','KCX':'LYS','B3E':'GLU','B3D':'ASP',
          'HZP':'PRO','CSX':'CYS','BAL':'ALA','HIC':'HIS',
          'DBZ':'ALA','DCY':'CYS','DVA':'VAL','NLE':'LEU',
          'SMC':'CYS','AGM':'ARG','B3A':'ALA','DAS':'ASP',
          'DLY':'LYS','DSN':'SER','DTH':'THR','GL3':'GLY',
          'HY3':'PRO','LLP':'LYS','MGN':'GLN','MHS':'HIS',
          'TRQ':'TRP','B3Y':'TYR','PHI':'PHE','PTR':'TYR',
          'TYS':'TYR','IAS':'ASP','GPL':'LYS','KYN':'TRP',
          'CSD':'CYS','SEC':'CYS'}

def pdb_to_string(pdb_file):
  modres = {**MODRES}
  lines = []
  for line in open(pdb_file,"rb"):
    line = line.decode("utf-8","ignore").rstrip()
    if line[:6] == "MODRES":
      k = line[12:15]
      v = line[24:27]
      if k not in modres and v in residue_constants.restype_3to1:
        modres[k] = v
    if line[:6] == "HETATM":
      k = line[17:20]
      if k in modres:
        line = "ATOM  "+line[6:17]+modres[k]+line[20:]
    if line[:4] == "ATOM":
      lines.append(line)
  return "\n".join(lines)
