import jax
import jax.numpy as jnp
import haiku as hk
import numpy as np

from .utils import cat_neighbors_nodes, get_ar_mask

class mpnn_score:
  def score(self, I):
    """ Graph-conditioned sequence model """
    
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

    if I.get("S",None) is None:
      ##########################################
      # unconditional_probs
      ##########################################      
      h_EXV_encoder_fw = h_EXV_encoder
      for layer in self.decoder_layers:
        h_V = layer(h_V, h_EXV_encoder_fw, I["mask"])
    else:
      ##########################################
      # conditional_probs
      ##########################################

      # Concatenate sequence embeddings for autoregressive decoder
      h_S = self.W_s(I["S"])
      h_ES = cat_neighbors_nodes(h_S, h_E, E_idx)

      if "ar_mask" not in I:
        if "decoding_order" not in I:
          key, sub_key = jax.random.split(key)
          randn = jax.random.uniform(sub_key, (I["X"].shape[0],))
          randn = jnp.where(I["mask"],randn,randn+1)
          if "fix_pos" in I: randn = randn.at[I["fix_pos"]].add(-1)
          I["decoding_order"] = randn.argsort()
        
        # make an autogressive mask
        I["ar_mask"] = get_ar_mask(I["decoding_order"])
              
      mask_attend = jnp.take_along_axis(I["ar_mask"], E_idx, 1)
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
    decoding_order = I.get("decoding_order",None)
    S = I.get("S",None)
    return {"logits": logits, "decoding_order":decoding_order, "S":S}
