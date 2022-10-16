import jax
import jax.numpy as jnp
import numpy as np
import re
import copy
import random
import os
import joblib
from tqdm import tqdm

from .modules import RunModel
from .utils import parse_PDB, StructureDatasetPDB, tied_featurize, _S_to_seq
from colabdesign.shared.prng import SafeKey
from colabdesign.mpnn.jax_weights import __file__ as mpnn_path

class MPNN_wrapper:  
  def __init__(self,
         model_name="v_48_020", verbose=False):
    self.model_name = model_name

    backbone_noise = 0.00  # Standard deviation of Gaussian noise to add to backbone atoms
    hidden_dim = 128
    num_layers = 3 

    path = os.path.join(os.path.dirname(mpnn_path), f'{model_name}.pkl')    
    checkpoint = joblib.load(path)
    params = jax.tree_util.tree_map(jnp.array, checkpoint['model_state_dict'])

    if verbose:
      print('Number of edges:', checkpoint['num_edges'])
      noise_level_print = checkpoint['noise_level']
      print(f'Training noise level: {noise_level_print}A')

    config = {'num_letters': 21,
          'node_features': hidden_dim,
          'edge_features': hidden_dim,
          'hidden_dim': hidden_dim,
          'num_encoder_layers': num_layers,
          'num_decoder_layers': num_layers,
          'augment_eps': backbone_noise,
          'k_neighbors': checkpoint['num_edges'],
          'dropout': 0.0
         }

    model = RunModel(config)
    model.params = params
    self.model = model

    self.alphabet = 'ACDEFGHIKLMNPQRSTVWYX'
    self.max_length = 20000

    seed = random.randint(0,2147483647)
    seed = jax.random.PRNGKey(seed)
    self.safe_key = SafeKey(seed)
  
  def prep_inputs(self, pdb_path,
          target_chain, fixed_chain=None,
          ishomomer=False, omit_AAs='X'):
    """generate input for score and sampling function

    Args:
      pdb_path (str): the path of the pdb file
      target_chain (str): chain ID of the protein sequence
      fixed_chain (str, optional): chain ID of the protein sequence that should be fixed. Defaults to None.
      ishomomer (bool, optional): for tie sampling. Defaults to False.
      omit_AAs (str, optional): aas should not be generated in sampling. Defaults to 'X'.

    Returns:
      dict: input dictionary
    """    
    # initialize some var
    fixed_positions_dict = None
    pssm_dict = None
    omit_AA_dict = None
    bias_by_res_dict = None
    bias_AAs_np = np.zeros(len(self.alphabet))
    pssm_threshold = 0.0
    pssm_multi = 0.0
    pssm_log_odds_flag = 0
    pssm_bias_flag = 0

    # fixed chain
    if fixed_chain is None:
      fixed_chain = ''
      fixed_chain_list = []
    else:
      fixed_chain_list = re.sub("[^A-Za-z]+",",", fixed_chain).split(",")

    # design chains
    if target_chain == '':
      designed_chain_list = []
    else:
      designed_chain_list = re.sub("[^A-Za-z]+",",", target_chain).split(",")

    #chain list
    chain_list = list(set(designed_chain_list + fixed_chain_list))

    # omit AAs
    omit_AAs_list = omit_AAs
    omit_AAs_np = np.array([AA in omit_AAs_list for AA in self.alphabet]).astype(np.float32)
  
    # prepare input
    pdb_dict_list = parse_PDB(pdb_path, input_chain_list=chain_list)
    dataset_valid = StructureDatasetPDB(pdb_dict_list, truncate=None, max_length=self.max_length)

    chain_id_dict = {}
    chain_id_dict[pdb_dict_list[0]['name']]= (designed_chain_list, fixed_chain_list)

    if ishomomer:
      # haven't tested
      tied_positions_dict = self.make_tied_positions_for_homomers(pdb_dict_list)
    else:
      tied_positions_dict = None

    return {'dataset_valid': dataset_valid,
        'chain_id_dict': chain_id_dict,
        'fixed_positions_dict': fixed_positions_dict,
        'omit_AA_dict': omit_AA_dict,
        'tied_positions_dict': tied_positions_dict,
        'pssm_dict': pssm_dict,
        'bias_by_res_dict': bias_by_res_dict,
        'pssm_threshold': pssm_threshold,
        'omit_AAs_np': omit_AAs_np,
        'bias_AAs_np': bias_AAs_np,
        'pssm_multi': pssm_multi,
        'pssm_log_odds_flag': pssm_log_odds_flag,
        'pssm_bias_flag': pssm_bias_flag,
         }
  
  def score(self, inputs, seq=None, order=None, key=None, unconditional=False):
    """get the output of MPNN

    Args:
      inputs (dict): output of the prep_input function
      seq (str, optional): the input sequence.
                 If not provided, the original sequence will be used.
                 Defaults to None.
      order (array, optional): the decoding order.
                   If not provided, the decoding order is random.
                   Defaults to None.
      key (jax.random.PRNGkey, optional): the random seed. Defaults to None.

    Returns:
      logits
      log_probs
    """    
    protein = inputs['dataset_valid'][0]
    batch_clones = [copy.deepcopy(protein)]
    (X, S, mask, lengths, chain_M, chain_idx, chain_list_list,
     visible_list_list, masked_list_list, masked_chain_length_list_list,
     chain_M_pos, omit_AA_mask, residue_idx, dihedral_mask,
     tied_pos_list_of_lists_list, pssm_coef, pssm_bias,
     pssm_log_odds_all, bias_by_res_all, tied_beta) = tied_featurize(batch_clones,
                                     inputs['chain_id_dict'], inputs['fixed_positions_dict'],
                                     inputs['omit_AA_dict'], inputs['tied_positions_dict'],
                                     inputs['pssm_dict'], inputs['bias_by_res_dict'])
    score_input = {'X': X,
                 'S': S,
                 'mask': mask,
                 'chain_M': chain_M * chain_M_pos,
                 'residue_idx': residue_idx,
                 'chain_idx': chain_idx}

    if unconditional:
      score_input["S"] = None
    else:
      if seq is not None:
        S = np.asarray([self.alphabet.index(a) for a in seq], dtype=np.int32)
        S = S[None, :]
        score_input['S'] = jnp.array(S)

      if order is None:
        if key is not None:
          self.safe_key = SafeKey(key)
        self.safe_key, used_key = self.safe_key.split()
        order = jax.random.normal(used_key.get(), (chain_M.shape[1],))
      score_input['randn'] = jnp.expand_dims(order, 0)
       
    self.safe_key, used_key = self.safe_key.split()
    return self.model.score(self.model.params, used_key.get(), score_input)
   
  def sampling(self, inputs,
         sample_num, batch_size,
         sampling_temp=0.1, order=None, key=None):
    """sample sequences from the given protein structure

    Args:
      inputs (dict): output of the prep_input function
      sample_num (int): number of sequences you want to generate
      batch_size (int): size of one batch
      sampling_temp (float, optional): sampling temperature. Defaults to 0.1.
      order (array, optional): the sampling order.
                   If not provided, the order is random.
                   Defaults to None.
      key (jax.random.PRNGkey, optional): the random seed. Defaults to None.

    Returns:
      seq_gen (list): generated sequence
    """
    NUM_BATCHES = sample_num//batch_size
    BATCH_COPIES = batch_size
    if key is not None:
      self.safe_key = SafeKey(key)

    protein = inputs['dataset_valid'][0]
    batch_clones = [copy.deepcopy(protein) for i in range(BATCH_COPIES)]
    (X, S, mask, lengths, chain_M, chain_idx, chain_list_list,
     visible_list_list, masked_list_list, masked_chain_length_list_list,
     chain_M_pos, omit_AA_mask, residue_idx, dihedral_mask,
     tied_pos_list_of_lists_list, pssm_coef, pssm_bias,
     pssm_log_odds_all, bias_by_res_all, tied_beta) = tied_featurize(batch_clones,
                                     inputs['chain_id_dict'], inputs['fixed_positions_dict'],
                                     inputs['omit_AA_dict'], inputs['tied_positions_dict'],
                                     inputs['pssm_dict'], inputs['bias_by_res_dict'])
    pssm_log_odds_mask = jax.lax.convert_element_type((pssm_log_odds_all > inputs['pssm_threshold']),
                              jnp.float32)  # 1.0 for true, 0.0 for false

    if order is None:
      self.safe_key, used_key = self.safe_key.split()
      order = jax.random.normal(used_key.get(), (chain_M.shape[1],))
    randn_1 = jnp.expand_dims(order, 0)

    # sample input
    sample_input = {'X': X,
            'randn': randn_1,
            'S_true': S,
            'chain_mask': chain_M,
            'chain_idx': chain_idx,
            'residue_idx': residue_idx,
            'mask': mask,
            'temperature': sampling_temp,
            'omit_AAs_np': inputs['omit_AAs_np'],
            'bias_AAs_np': inputs['bias_AAs_np'],
            'chain_M_pos': chain_M_pos,
            'omit_AA_mask': omit_AA_mask,
            'pssm_coef': pssm_coef,
            'pssm_bias': pssm_bias,
            'pssm_multi': inputs['pssm_multi'],
            'pssm_log_odds_flag': bool(inputs['pssm_log_odds_flag']),
            'pssm_log_odds_mask': pssm_log_odds_mask,
            'pssm_bias_flag': bool(inputs['pssm_bias_flag']),
            'bias_by_res': bias_by_res_all
            }
    seq_gen = []
    for _ in tqdm(range(NUM_BATCHES)):
      self.safe_key, used_key = self.safe_key.split()
      sample_input.update({'key': used_key.get()})

      self.safe_key, used_key = self.safe_key.split()
      if inputs['tied_positions_dict'] is None:
        sample_dict = self.model.sample(self.model.params, used_key.get(), sample_input)
      else:
        sample_input.update({'tied_pos': tied_pos_list_of_lists_list[0],
                   'tied_beta': tied_beta,
                   'bias_by_res': bias_by_res_all,
                  })
        sample_dict = self.model.tied_sample(self.model.params, used_key.get(), sample_input)
      S_sample = sample_dict["S"]
      for b_ix in range(BATCH_COPIES):
        masked_chain_length_list = masked_chain_length_list_list[b_ix]
        masked_list = masked_list_list[b_ix]
        seq = _S_to_seq(S_sample[b_ix], chain_M[b_ix])

        start = 0
        end = 0
        list_of_AAs = []
        for mask_l in masked_chain_length_list:
          end += mask_l
          list_of_AAs.append(seq[start:end])
          start = end

        seq = "".join(list(np.array(list_of_AAs)[np.argsort(masked_list)]))
        l0 = 0
        for mc_length in list(np.array(masked_chain_length_list)[np.argsort(masked_list)])[:-1]:
          l0 += mc_length
          seq = seq[:l0] + '/' + seq[l0:]
          l0 += 1
        seq_gen.append(seq)
    return seq_gen

  @staticmethod
  def make_tied_positions_for_homomers(pdb_dict_list):
    my_dict = {}
    for result in pdb_dict_list:
      all_chain_list = sorted([item[-1:] for item in list(result) if item[:9]=='seq_chain'])  # A, B, C, ...
      tied_positions_list = []
      chain_length = len(result[f"seq_chain_{all_chain_list[0]}"])
      for i in range(1,chain_length+1):
        temp_dict = {}
        for j, chain in enumerate(all_chain_list):
          temp_dict[chain] = [i] #needs to be a list
        tied_positions_list.append(temp_dict)
      my_dict[result['name']] = tied_positions_list
    return my_dict
