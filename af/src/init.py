import tensorflow as tf
tf.config.set_visible_devices([], 'GPU')
import jax
import jax.numpy as jnp
import numpy as np

from alphafold.data import pipeline, prep_inputs
from alphafold.common import protein, residue_constants
from alphafold.model import all_atom
from alphafold.model.tf import shape_placeholders

from af.src.misc import _np_get_cb
ORDER_RESTYPE = {v: k for k, v in residue_constants.restype_order.items()}

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
      batch.update(all_atom.atom37_to_frames(**batch))

      seq = "".join([ORDER_RESTYPE[a] for a in batch["aatype"]])
      template_aatype = residue_constants.sequence_to_onehot(seq, residue_constants.HHBLITS_AA_TO_ID)
      template_features = {"template_aatype":template_aatype,
                           "template_all_atom_masks":batch["all_atom_mask"],
                           "template_all_atom_positions":batch["all_atom_positions"]}
      
      residue_index = protein_obj.residue_index[has_ca] + last      
      last = residue_index[-1] + 50
      o.append({"batch":batch,
                "template_features":template_features,
                "residue_index": residue_index})
      
      # original residue and chain index
      residue_idx.append(protein_obj.residue_index[has_ca])
      chain_idx.append([chain] * len(residue_idx[-1]))

    o = jax.tree_util.tree_map(lambda *x:np.concatenate(x,0),*o)
    o["template_features"] = {"template_domain_names":np.asarray(["None"]),
                              **jax.tree_map(lambda x:x[None],o["template_features"])}

    # save original residue and chain index
    o["idx"] = {"residue":np.concatenate(residue_idx),
                "chain":np.concatenate(chain_idx)}
    return o
  
  def _prep_pos(self, pos, residue=None, chain=None):
    if residue is None: residue = self._inputs["residue_index"][0]
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
  
  # prep functions specific to protocol
  def _prep_binder(self, pdb_filename, chain="A",
                   binder_len=50, binder_chain=None, use_binder_template=False,
                   hotspot=None, **kwargs):

    '''prep inputs for binder design'''
    
    self._redesign = binder_chain is not None

    # get pdb info
    chains = f"{chain},{binder_chain}" if self._redesign else chain
    pdb = self._prep_pdb(pdb_filename, chain=chains)    
    if self._redesign:
      target_len = sum([(pdb["idx"]["chain"] == c).sum() for c in chain.split(",")])
      binder_len = sum([(pdb["idx"]["chain"] == c).sum() for c in binder_chain.split(",")])
      self._inputs = self._prep_features(target_len + binder_len, pdb["template_features"])
    else:
      target_len = pdb["residue_index"].shape[0]
      self._inputs = self._prep_features(target_len, pdb["template_features"])
      
    self._inputs["residue_index"][...,:] = pdb["residue_index"]

    self._copies = 1
    self._repeat = False

    # gather hotspot info
    if hotspot is None:
      self._hotspot = None
    else:
      self._hotspot = self._prep_pos(hotspot, **pdb["idx"])

    if self._redesign:
      self._inputs["template_aatype"][...,target_len:] = 21
      self._inputs["template_all_atom_masks"][...,target_len:,5:] = 0
      if not use_binder_template:
        self._inputs["template_all_atom_masks"][...,target_len:,:] = 0
        self._inputs["template_pseudo_beta_mask"][...,target_len:] = 0      
      
      self._batch = pdb["batch"]
      self._wt_aatype = self._batch["aatype"][target_len:]
      self._default_weights.update({"dgram_cce":1.0, "fape":0.0, "rmsd":0.0,
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
      self._default_weights.update({"con":0.5, "i_pae":0.01, "i_con":0.5})

    self._target_len = target_len
    self._binder_len = self._len = binder_len

    self.restart(set_defaults=True, **kwargs)

  def _prep_fixbb(self, pdb_filename, chain=None, copies=1, homooligomer=False, repeat=False, **kwargs):
    '''prep inputs for fixed backbone design'''
    pdb = self._prep_pdb(pdb_filename, chain=chain)
    self._batch = pdb["batch"]
    self._wt_aatype = self._batch["aatype"]
    self._len = pdb["residue_index"].shape[0]
    temp_features = {'template_aatype': np.zeros([1,self._len,22]),
                     'template_all_atom_masks': np.zeros([1,self._len,37]),
                     'template_all_atom_positions': np.zeros([1,self._len,37,3]),
                     'template_domain_names': np.asarray(["None"])}
    self._inputs = self._prep_features(self._len, temp_features)
    self._copies = copies
    self._repeat = repeat
    
    # set weights
    self._default_weights.update({"dgram_cce":1.0, "fape":0.0, "rmsd":0.0, "con":0.0})
    self._default_opt.update({"6D_weights":{}})

    # update residue index from pdb
    if copies > 1:
      if repeat:
        self._len = self._len // copies
      else:
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
    
  def _prep_hallucination(self, length=100, copies=1, repeat=False, **kwargs):
    '''prep inputs for hallucination'''
    self._len = length
    self._inputs = self._prep_features(length * copies)
    self._copies = copies
    self._repeat = repeat
    
    # set weights
    self._default_weights.update({"con":1.0})
    if copies > 1:
      if repeat:
        offset = 1
      else:
        offset = 50
        self._default_weights.update({"i_pae":0.01, "i_con":0.1})
      self._inputs["residue_index"] = self.repeat_idx(np.arange(length), copies, offset=offset)[None]

    self.restart(set_defaults=True, **kwargs)

  def _prep_partial(self, pdb_filename, chain=None, pos=None, length=None,
                    fix_seq=True, sidechain=False, **kwargs):
    '''prep input for partial hallucination'''
    self.args.update({"sidechain":sidechain, "fix_seq":fix_seq})
    self._copies = 1
    self._repeat = False    
    
    # get [pos]itions of interests
    pdb = self._prep_pdb(pdb_filename, chain=chain)
    if pos is None:
      pos = np.arange(pdb["residue_index"].shape[0])
    else:
      pos = self._prep_pos(pos, **pdb["idx"])
    
    self._default_opt["pos"] = pos
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
    else:
      self._inputs = self._prep_features(self._len)
    
    weights = {"dgram_cce":1.0,"con":1.0, "fape":0.0, "rmsd":0.0,"6D":0.0}
    if sidechain: weights.update({"sc_fape":0.0, "sc_rmsd":0.0})
    self._default_weights.update(weights)

    self.restart(set_defaults=True, **kwargs)

#######################
# utils
#######################
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
  lines = []
  for line in open(pdb_file,"r"):
    if line[:6] == "HETATM" and line[17:20] in MODRES:
      line = "ATOM  "+line[6:17]+MODRES[line[17:20]]+line[20:]
    if line[:4] == "ATOM":
      lines.append(line)
  return "".join(lines)
