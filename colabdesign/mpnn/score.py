import jax
import jax.numpy as jnp
import haiku as hk
import numpy as np

from .utils import cat_neighbors_nodes, get_ar_mask

class mpnn_score:
  def score(self, I):
    """
    I = {
         [[required]]
         'X' = (L,4,3) 
         'mask' = (L,)
         'residue_index' = (L,)
         'chain_idx' = (L,)
         
         [[optional]]
         'S' = (L,21)
         'decoding_order' = (L,)
         'ar_mask' = (L,L)
        }
    """
    
    key = hk.next_rng_key()
    # Prepare node and edge embeddings
    E, E_idx = self.features(I)
    h_V = jnp.zeros((E.shape[0], E.shape[-1]))
    h_E = self.W_e(E)
    
    # Encoder is unmasked self-attention
    mask_attend = jnp.take_along_axis(I["mask"][:,None] * I["mask"][None,:], E_idx, 1)
    
    for layer in self.encoder_layers:
      h_V, h_E = layer(h_V, h_E, E_idx, I["mask"], mask_attend)

    # Build encoder embeddings
    h_EX_encoder = cat_neighbors_nodes(jnp.zeros_like(h_V), h_E, E_idx)
    h_EXV_encoder = cat_neighbors_nodes(h_V, h_EX_encoder, E_idx)

    if "S" not in I:
      ##########################################
      # unconditional_probs
      ##########################################      
      h_EXV_encoder_fw = h_EXV_encoder
      for layer in self.decoder_layers:
        h_V = layer(h_V, h_EXV_encoder_fw, I["mask"])
      decoding_order = None
    else:
      ##########################################
      # conditional_probs
      ##########################################

      # Concatenate sequence embeddings for autoregressive decoder
      h_S = self.W_s(I["S"])
      h_ES = cat_neighbors_nodes(h_S, h_E, E_idx)

      # get autoregressive mask
      if "ar_mask" in I:
        decoding_order = None
        ar_mask = I["ar_mask"]
      else:
        decoding_order = I["decoding_order"]
        ar_mask = get_ar_mask(decoding_order)
              
      mask_attend = jnp.take_along_axis(ar_mask, E_idx, 1)
      mask_1D = I["mask"][:,None]
      mask_bw = mask_1D * mask_attend
      mask_fw = mask_1D * (1 - mask_attend)

      h_EXV_encoder_fw = mask_fw[...,None] * h_EXV_encoder
      for layer in self.decoder_layers:
        # Masked positions attend to encoder information, unmasked see. 
        h_ESV = cat_neighbors_nodes(h_V, h_ES, E_idx)
        h_ESV = mask_bw[...,None] * h_ESV + h_EXV_encoder_fw
        h_V = layer(h_V, h_ESV, I["mask"])
    
    logits = self.W_out(h_V)
    S = I.get("S",None)
    return {"logits": logits, "decoding_order":decoding_order, "S":S}