# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
confidence.py, but with jax instead of numpy
converted by @konstin
"""

import jax
from jax import jit

def compute_plddt_jax(logits):
    """Port of confidence.compute_plddt to jax

    Computes per-residue pLDDT from logits.

    Args:
      logits: [num_res, num_bins] output from the PredictedLDDTHead.

    Returns:
      plddt: [num_res] per-residue pLDDT.
    """
    num_bins = logits.shape[-1]
    bin_width = 1.0 / num_bins
    bin_centers = jax.numpy.arange(start=0.5 * bin_width, stop=1.0, step=bin_width)
    probs = jax.nn.softmax(logits, axis=-1)
    predicted_lddt_ca = jax.numpy.sum(probs * bin_centers[None, :], axis=-1)
    return predicted_lddt_ca * 100

def _calculate_bin_centers(breaks):
    """Gets the bin centers from the bin edges.

    Args:
      breaks: [num_bins - 1] the error bin edges.

    Returns:
      bin_centers: [num_bins] the error bin centers.
    """
    step = (breaks[1] - breaks[0])

    # Add half-step to get the center
    bin_centers = breaks + step / 2
    # Add a catch-all bin at the end.
    bin_centers = jax.numpy.concatenate([bin_centers, jax.numpy.asarray([bin_centers[-1] + step])],
                                        axis=0)
    return bin_centers

def predicted_tm_score_jax(
    logits,
    breaks,
    residue_weights=None,
    asym_id=None,
    interface=False):
    """Computes predicted TM alignment or predicted interface TM alignment score.

    Args:
      logits: [num_res, num_res, num_bins] the logits output from
        PredictedAlignedErrorHead.
      breaks: [num_bins] the error bins.
      residue_weights: [num_res] the per residue weights to use for the
        expectation.
      asym_id: [num_res] the asymmetric unit ID - the chain ID. Only needed for
        ipTM calculation, i.e. when interface=True.
      interface: If True, interface predicted TM score is computed.

    Returns:
      ptm_score: The predicted TM alignment or the predicted iTM score.
    """

    # residue_weights has to be in [0, 1], but can be floating-point, i.e. the
    # exp. resolved head's probability.
    if residue_weights is None:
        residue_weights = jax.numpy.ones(logits.shape[0])

    bin_centers = _calculate_bin_centers(breaks)

    num_res = logits.shape[0] # jax.numpy.sum(residue_weights)
    # Clip num_res to avoid negative/undefined d0.
    clipped_num_res = jax.numpy.max(jax.numpy.asarray([num_res, 19]))

    # Compute d_0(num_res) as defined by TM-score, eqn. (5) in Yang & Skolnick
    # "Scoring function for automated assessment of protein structure template
    # quality", 2004: http://zhanglab.ccmb.med.umich.edu/papers/2004_3.pdf
    d0 = 1.24 * (clipped_num_res - 15) ** (1. / 3) - 1.8

    # Convert logits to probs.
    probs = jax.nn.softmax(logits, axis=-1)

    # TM-Score term for every bin.
    tm_per_bin = 1. / (1 + jax.numpy.square(bin_centers) / jax.numpy.square(d0))
    # E_distances tm(distance).
    predicted_tm_term = jax.numpy.sum(probs * tm_per_bin, axis=-1)

    pair_mask = jax.numpy.ones(shape=(num_res, num_res), dtype=bool)
    if interface:
        pair_mask *= asym_id[:, None] != asym_id[None, :]

    predicted_tm_term *= pair_mask

    pair_residue_weights = pair_mask * (
        residue_weights[None, :] * residue_weights[:, None])
    normed_residue_mask = pair_residue_weights / (1e-8 + pair_residue_weights.sum(-1,keepdims=True))
    per_alignment = jax.numpy.sum(predicted_tm_term * normed_residue_mask, axis=-1)
    return jax.numpy.asarray(per_alignment[(per_alignment * residue_weights).argmax()])

def get_confidence_metrics(
    prediction_result,
    multimer_mode: bool):
    """Post processes prediction_result to get confidence metrics."""
    confidence_metrics = {}
    confidence_metrics['plddt'] = compute_plddt_jax(
        prediction_result['predicted_lddt']['logits'])
    if 'predicted_aligned_error' in prediction_result:
        confidence_metrics.update(compute_predicted_aligned_error(
            logits=prediction_result['predicted_aligned_error']['logits'],
            breaks=prediction_result['predicted_aligned_error']['breaks']))
        confidence_metrics['ptm'] = predicted_tm_score_jax(
            logits=prediction_result['predicted_aligned_error']['logits'],
            breaks=prediction_result['predicted_aligned_error']['breaks'],
            asym_id=None)
        if multimer_mode:
            # Compute the ipTM only for the multimer model.
            confidence_metrics['iptm'] = predicted_tm_score_jax(
                logits=prediction_result['predicted_aligned_error']['logits'],
                breaks=prediction_result['predicted_aligned_error']['breaks'],
                asym_id=prediction_result['predicted_aligned_error']['asym_id'],
                interface=True)
            confidence_metrics['ranking_confidence'] = (
                0.8 * confidence_metrics['iptm'] + 0.2 * confidence_metrics['ptm'])

    if not multimer_mode:
        # Monomer models use mean pLDDT for model ranking.
        confidence_metrics['ranking_confidence'] = jax.numpy.mean(
            confidence_metrics['plddt'])

    return confidence_metrics

def _calculate_expected_aligned_error(
    alignment_confidence_breaks,
    aligned_distance_error_probs):
  """Calculates expected aligned distance errors for every pair of residues.

  Args:
    alignment_confidence_breaks: [num_bins - 1] the error bin edges.
    aligned_distance_error_probs: [num_res, num_res, num_bins] the predicted
      probs for each error bin, for each pair of residues.

  Returns:
    predicted_aligned_error: [num_res, num_res] the expected aligned distance
      error for each pair of residues.
    max_predicted_aligned_error: The maximum predicted error possible.
  """
  bin_centers = _calculate_bin_centers(alignment_confidence_breaks)

  # Tuple of expected aligned distance error and max possible error.
  return (jax.numpy.sum(aligned_distance_error_probs * bin_centers, axis=-1),
          jax.numpy.asarray(bin_centers[-1]))

def compute_predicted_aligned_error(
    logits,
    breaks):
    """Computes aligned confidence metrics from logits.

    Args:
      logits: [num_res, num_res, num_bins] the logits output from
        PredictedAlignedErrorHead.
      breaks: [num_bins - 1] the error bin edges.

    Returns:
      aligned_confidence_probs: [num_res, num_res, num_bins] the predicted
        aligned error probabilities over bins for each residue pair.
      predicted_aligned_error: [num_res, num_res] the expected aligned distance
        error for each pair of residues.
      max_predicted_aligned_error: The maximum predicted error possible.
    """
    aligned_confidence_probs = jax.nn.softmax(
        logits,
        axis=-1)
    predicted_aligned_error, max_predicted_aligned_error = (
        _calculate_expected_aligned_error(
            alignment_confidence_breaks=breaks,
            aligned_distance_error_probs=aligned_confidence_probs))
    return {
        'aligned_confidence_probs': aligned_confidence_probs,
        'predicted_aligned_error': predicted_aligned_error,
        'max_predicted_aligned_error': max_predicted_aligned_error,
    }
