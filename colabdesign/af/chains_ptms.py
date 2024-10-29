import numpy as np
from colabdesign.af.alphafold.common import confidence
from colabdesign.af import loss
import string
from copy import deepcopy
import jax.numpy as jnp
import jax
import scipy

""" Functions to calculate interface metrics on actual interfaces"""


def get_chain_indices(chain_boundaries):
    """Returns a list of tuples indicating the start and end indices for each chain."""
    chain_starts_ends = []
    unique_chains = np.unique(chain_boundaries)
    for chain in unique_chains:
        positions = np.where(chain_boundaries == chain)[0]
        chain_starts_ends.append((positions[0], positions[-1]))
    return chain_starts_ends


def predicted_tm_score_modified(logits, breaks, residue_weights=None,
                                asym_id=None, use_jnp=False, pair_residue_weights=None):
    """Computes predicted TM alignment or predicted interface TM alignment score.

    Args:
      logits: [num_res, num_res, num_bins] the logits output from
        PredictedAlignedErrorHead.
      breaks: [num_bins] the error bins.
      residue_weights: [num_res] the per residue weights to use for the
        expectation.
      asym_id: [num_res] the asymmetric unit ID - the chain ID. Only needed for
        ipTM calculation.
      pair_residue_weights: [num_res, num_res] unnormalized weights for actifptm calculation

    Returns:
      ptm_score: The predicted TM alignment or the predicted iTM score.
    """
    if use_jnp:
        _np, _softmax = jnp, jax.nn.softmax
    else:
        _np, _softmax = np, scipy.special.softmax

    # residue_weights has to be in [0, 1], but can be floating-point, i.e. the
    # exp. resolved head's probability.
    if residue_weights is None:
        residue_weights = _np.ones(logits.shape[0])
    bin_centers = confidence._calculate_bin_centers(breaks, use_jnp=use_jnp)
    num_res = residue_weights.shape[0]

    # Clip num_res to avoid negative/undefined d0.
    clipped_num_res = _np.maximum(num_res, 19)

    # Compute d_0(num_res) as defined by TM-score, eqn. (5) in Yang & Skolnick
    # "Scoring function for automated assessment of protein structure template
    # quality", 2004: http://zhanglab.ccmb.med.umich.edu/papers/2004_3.pdf
    d0 = 1.24 * (clipped_num_res - 15) ** (1. / 3) - 1.8

    # Convert logits to probs.
    probs = _softmax(logits, axis=-1)

    # TM-Score term for every bin.
    tm_per_bin = 1. / (1 + _np.square(bin_centers) / _np.square(d0))
    import seaborn as sns
    import matplotlib.pyplot as plt
    # E_distances tm(distance).
    predicted_tm_term = (probs * tm_per_bin).sum(-1)

    if asym_id is None:
        pair_mask = _np.full((num_res, num_res), True)
    else:
        pair_mask = asym_id[:, None] != asym_id[None, :]

    predicted_tm_term *= pair_mask

    # If pair_residue_weights is provided (e.g. for if_ptm with contact probabilities),
    # it should not be overwritten
    if pair_residue_weights is None:
        pair_residue_weights = pair_mask * (residue_weights[None, :] * residue_weights[:, None])
    normed_residue_mask = pair_residue_weights / (1e-8 + pair_residue_weights.sum(-1, keepdims=True))

    per_alignment = (predicted_tm_term * normed_residue_mask).sum(-1)

    return (per_alignment * residue_weights).max()


def get_ptm_modified(inputs, outputs, interface=False):
    """ This function is the same as in the original AF2, just calls a modified TM-score calculation function."""

    pae = {"residue_weights": inputs["seq_mask"], **outputs["predicted_aligned_error"]}

    if interface:
        if "asym_id" not in pae:
            pae["asym_id"] = inputs["asym_id"]
    else:
        if "asym_id" in pae:
            pae.pop("asym_id")

    return predicted_tm_score_modified(**pae, use_jnp=True)


def get_actifptm_probs(af, cmap, start_i, end_i, start_j, end_j):
    """
    This function calculates the interface PTM score for a given interface, taking into account the contact probabilities, not binary contacts.

    Args:
        af: AlphaFold object
        cmap: Contact map
        start_i, end_i: Start and end indices of the first chain
        start_j, end_j: Start and end indices of the second chain

    Returns:
        actifptm: Interface pTM score
        ipTM: Interchain pTM score
    """

    total_length = len(af._inputs['asym_id'])
    outputs = deepcopy(af.aux['debug']['outputs'])
    inputs_actifptm = {}

    # Create a new matrix, which contains only the contacts between the two chains
    cmap_copy = np.zeros((total_length, total_length))
    cmap_copy[start_i:end_i + 1, start_j:end_j + 1] = cmap[start_i:end_i + 1, start_j:end_j + 1]
    cmap_copy[start_j:end_j + 1, start_i:end_i + 1] = cmap[start_j:end_j + 1, start_i:end_i + 1]

    # Initialize new input dictionary
    inputs_actifptm['seq_mask'] = np.full(total_length, 0, dtype=float)
    inputs_actifptm['asym_id'] = af._inputs['asym_id']
    # Update seq_mask for these positions to True within inputs
    inputs_actifptm['seq_mask'][np.concatenate((np.arange(start_i, end_i + 1),
                                             np.arange(start_j, end_j + 1)))] = 1

    # Call get_ptm with updated inputs and outputs
    pae = {"residue_weights": inputs_actifptm["seq_mask"],
           **outputs["predicted_aligned_error"],
           "asym_id": inputs_actifptm["asym_id"]}
    actifptm_probs = round(float(predicted_tm_score_modified(**pae, use_jnp=True,
                                                          pair_residue_weights=cmap_copy)), 3)

    return actifptm_probs


def get_actifptm_contacts(af, cmap, start_i, end_i, start_j, end_j):
    """
    This function calculates the interface PTM score for a given interface, taking into account binary contacts.
    In case of no confident contacts, the interface PTM score is set to 0.

    Args:
        af: AlphaFold object
        cmap: Contact map
        start_i, end_i: Start and end indices of the first chain
        start_j, end_j: Start and end indices of the second chain

    Returns:
        actifptm: Interface pTM score
        ipTM: Interchain pTM score
    """

    # Prepare a dictionary to collect the inputs for calculation
    inputs_actifptm = {}
    contacts = np.where(cmap[start_i:end_i + 1, start_j:end_j + 1] >= 0.6)
    total_length = len(af._inputs['asym_id'])
    outputs = deepcopy(af.aux['debug']['outputs'])

    if contacts[0].size > 0:  # If there are contacts
        # Convert local chain positions back to global positions using JAX
        global_i_positions = contacts[0] + start_i
        global_j_positions = contacts[1] + start_j
        global_positions = list(set(np.concatenate((global_i_positions, global_j_positions))))
        global_positions = np.array(global_positions, dtype=int)
        global_positions.sort()

        # Initialize new input dictionary
        inputs_actifptm['seq_mask'] = np.full(total_length, 0, dtype=float)
        inputs_actifptm['asym_id'] = af._inputs['asym_id']
        # Update seq_mask for these positions to True within inputs
        inputs_actifptm['seq_mask'][global_positions] = 1
        # Call get_ptm with updated inputs and outputs
        actifptm = round(float(get_ptm_modified(inputs_actifptm, outputs, interface=True)), 3)
    else:
        actifptm = 0

    return actifptm


def get_pairwise_iptm(af, start_i, end_i, start_j, end_j):
    """This will calculate ipTM as usual, just between given chains"""

    input_pairwise_iptm = {}

    # Prepare inputs
    outputs = deepcopy(af.aux['debug']['outputs'])
    input_pairwise_iptm['seq_mask'] = np.full(len(af._inputs['asym_id']), 0, dtype=float)
    input_pairwise_iptm['asym_id'] = af._inputs['asym_id']
    input_pairwise_iptm['seq_mask'][np.concatenate((np.arange(start_i, end_i + 1),
                                                    np.arange(start_j, end_j + 1)))] = 1

    return round(float(get_ptm_modified(input_pairwise_iptm, outputs, interface=True)), 3)


def get_per_chain_ptm(af, cmap, start, end):
    """
    Calculates the chain PTM score for a specified interface region, using contact probabilities.

    Args:
        af: AlphaFold object
        cmap: Contact map
        start: Start index of the interface region
        end: End index of the interface region

    Returns:
        cpTM: Chain pTM score
    """
    # Extract only the relevant subset of the contact map
    cmap_subset = cmap[start:end + 1, start:end + 1]

    # Extract only the relevant subset of the logits map
    pae_copy = deepcopy(af.aux['debug']['outputs'])["predicted_aligned_error"]
    pae_copy['logits'] = pae_copy['logits'][start:end + 1, start:end + 1, :]

    # Prepare inputs for the modified predicted_tm_score function
    pae_copy['residue_weights'] = np.ones(end - start + 1, dtype=float)
    del (pae_copy['asym_id'])

    # Calculate and return chain PTM score
    cptm_probs = round(float(predicted_tm_score_modified(**pae_copy, use_jnp=True)), 3)

    return cptm_probs


def get_pairwise_interface_metrics(af, calculate_interface=False):
    """
    This function iterates over all pairs of chains and calculates the interface and interchain PTM score for each pair.

    Args:
        af: AlphaFold object
        calculate_interface: If True, calculate interface pTM score based on contact probabilities. If False, calculate interface pTM score based on binary contacts.

    Returns:
        pairwise_actifptm: Dictionary containing pairwise interface pTM scores
        pairwise_iptm: Dictionary containing pairwise interchain pTM scores
    """

    # Prepare a dictionary to collect results
    pairwise_actifptm = {}
    pairwise_iptm = {}
    per_chain_ptm = {}
    chain_starts_ends = get_chain_indices(af._inputs['asym_id'])

    # Generate chain labels (A, B, C, ...)
    chain_labels = list(string.ascii_uppercase)

    cmap = loss.get_contact_map(af.aux['debug']['outputs'], 8)  # Define interface with 8A between Cb-s
    for i, (start_i, end_i) in enumerate(chain_starts_ends):
        chain_label_i = chain_labels[i % len(chain_labels)]  # Wrap around if more than 26 chains
        for j, (start_j, end_j) in enumerate(chain_starts_ends):
            chain_label_j = chain_labels[j % len(chain_labels)]  # Wrap around if more than 26 chains
            if i < j:  # Avoid self-comparison and duplicate comparisons
                key = f"{chain_label_i}-{chain_label_j}"

                if calculate_interface:
                    pairwise_actifptm[key] = get_actifptm_contacts(af, cmap, start_i, end_i, start_j, end_j)
                else:
                    pairwise_actifptm[key] = get_actifptm_probs(af, cmap, start_i, end_i, start_j, end_j)

                # Also add regular i_ptm (interchain), pairwise
                pairwise_iptm[key] = get_pairwise_iptm(af, start_i, end_i, start_j, end_j)

        # Also calculate pTM score for single chain
        per_chain_ptm[chain_label_i] = get_per_chain_ptm(af, cmap, start_i, end_i)

    return pairwise_actifptm, pairwise_iptm, per_chain_ptm
