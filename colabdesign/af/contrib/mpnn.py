import jax
import jax.numpy as jnp

from colabdesign.mpnn import mk_mpnn_model
from colabdesign.af.alphafold.common import residue_constants

def add_mpnn_loss(self, mode="ar_mask", entropy=False):
  '''add mpnn loss'''    
  
  def loss_fn(inputs, outputs, aux, params, key):

    # get structure
    atom_idx = tuple(residue_constants.atom_order[k] for k in ["N","CA","C","O"])    
    I = {"S":           inputs["aatype"],
         "mask":        inputs["seq_mask"],
         "residue_idx": inputs["residue_index"],
         "chain_idx":   inputs["asym_id"],
         "X":           outputs["structure_module"]["final_atom_positions"][:,atom_idx],
         "lengths":     self._lengths,
         "key":         key}
    
    # add autoregressive mask
    if mode == "ar_mask":
      L = sum(self._lengths)
      fix_pos = inputs["fix_pos"]
      fix_pos = jnp.logical_or(fix_pos[:,None],fix_pos[None,:])
      I["ar_mask"] = jnp.where(fix_pos,1-jnp.eye(L),0)
    else:
      I["fix_pos"] = inputs["fix_pos"]

    # adjust offset
    if "offset" in inputs: I["offset"] = inputs["offset"]

    # get logits
    aux["mpnn_logits"] = self._mpnn._score(**I)["logits"][:,:20]
    if entropy:
      # entropy
      q_probs = jax.nn.softmax(aux["mpnn_logits"])
      q_logprobs = jax.nn.log_softmax(aux["mpnn_logits"])
      loss = -(q_probs * q_logprobs).sum(-1)
    
    else:
      # kl_divergence
      q_logits = jax.lax.stop_gradient(aux["mpnn_logits"])      
      q_probs = jax.nn.softmax(q_logits)
      q_logprobs = jax.nn.log_softmax(q_logits)
      p_logprobs = jax.nn.log_softmax(params["seq"]*inputs["opt"]["alpha"])      
      loss_kl = (q_probs * (q_logprobs - p_logprobs)).sum(-1)

      # cross_entropy
      p = inputs["msa"][...,:20]
      p_hard = jax.nn.one_hot(p.argmax(-1),20)
      p = jax.lax.stop_gradient(p_hard - p) + p
      loss_cce = -(p.mean(0) * q_logprobs).sum(-1)
      
      # switch
      loss = jnp.where(inputs["opt"]["mpnn_kl"],loss_kl,loss_cce)
    
    mask = jnp.where(inputs["fix_pos"],False,inputs["seq_mask"])
    return {"mpnn":(loss * mask).sum() / (mask.sum() + 1e-8)}
  
  self._mpnn = mk_mpnn_model()  
  self._callbacks["model"]["loss"].append(loss_fn)
  self._opt["weights"]["mpnn"] = 1.0
  
  if not entropy:
    self._opt["mpnn_kl"] = False