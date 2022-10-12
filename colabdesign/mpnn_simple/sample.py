import jax
import jax.numpy as jnp
import numpy as np
import itertools

from colabdesign.shared.prng import SafeKey
from .utils import gather_nodes, cat_neighbors_nodes, scatter, get_ar_mask

class mpnn_sample:
  def sample(self, X, decoding_order,
             chain_idx, residue_idx,
             mask=None, temperature=1.0):

    key = SafeKey(hk.next_rng_key())

    # Prepare node and edge embeddings
    E, E_idx = self.features(X, mask, residue_idx, chain_idx)
    h_V = jnp.zeros((E.shape[0], E.shape[1], E.shape[-1]))
    h_E = self.W_e(E)

    # Encoder is unmasked self-attention
    mask_attend = gather_nodes(mask[...,None], E_idx)[...,0]
    mask_attend = mask[...,None] * mask_attend
    for layer in self.encoder_layers:
      h_V, h_E = layer(h_V, h_E, E_idx, mask, mask_attend)

    # Decoder uses masked self-attention
    ar_mask = get_ar_mask(decoding_order)

    mask_attend = jnp.take_along_axis(ar_mask, E_idx, 2)[...,None]
    mask_1D = mask.reshape([mask.shape[0], mask.shape[1], 1, 1])
    mask_bw = mask_1D * mask_attend
    mask_fw = mask_1D * (1. - mask_attend)

    N_batch, N_nodes = X.shape[0], X.shape[1]
    
    # initial values
    all_logits = jnp.zeros((N_batch, N_nodes, 21))
    h_S = jnp.zeros_like(h_V,)
    S = jnp.zeros((N_batch, N_nodes),float)    
    h_V_stack = [h_V] + [jnp.zeros_like(h_V) for _ in range(len(self.decoder_layers))]
    h_EX_encoder = cat_neighbors_nodes(jnp.zeros_like(h_S), h_E, E_idx)
    h_EXV_encoder = cat_neighbors_nodes(h_V, h_EX_encoder, E_idx)
    h_EXV_encoder_fw = mask_fw * h_EXV_encoder

    def step(i):
      t = decoding_order[:,i]  # [B]
      t1 = t[:,None]
      t2 = t[:,None,None]
      t3 = t[:,None,None,None]
      
      # Hidden layers
      t_ = jnp.tile(t2, [1,1,E_idx.shape[-1]])
      E_idx_t = jnp.take_along_axis(E_idx, t_, 1)
      t_ = jnp.tile(t3, [1,1,h_E.shape[-2], h_E.shape[-1]])
      h_E_t = jnp.take_along_axis(h_E, t_, 1)
      h_ES_t = cat_neighbors_nodes(h_S, h_E_t, E_idx_t)
      t_ = jnp.tile(t3, [1,1,h_EXV_encoder_fw.shape[-2], h_EXV_encoder_fw.shape[-1]])
      h_EXV_encoder_t = jnp.take_along_axis(h_EXV_encoder_fw, t_, 1)
      mask_t = jnp.take_along_axis(mask, t1, 1)

      for l, layer in enumerate(self.decoder_layers):
        h_V_stack_ = h_V_stack[l]
        # Updated relational features for future states
        h_ESV_decoder_t = cat_neighbors_nodes(h_V_stack_, h_ES_t, E_idx_t)
        t_ = jnp.tile(t2, [1,1,h_V_stack_.shape[-1]])
        h_V_t = jnp.take_along_axis(h_V_stack_, t_, 1)
        t_ = jnp.tile(t3, [1,1,mask_bw.shape[-2], mask_bw.shape[-1]])
        mask_bw_ = jnp.take_along_axis(mask_bw, t_,1)
        h_ESV_t = mask_bw_  * h_ESV_decoder_t + h_EXV_encoder_t

        # updates
        t_ = jnp.tile(t2 [1,1,h_V.shape[-1]])
        h_V_stack[l+1] = scatter(h_V_stack[l+1], 1, t_, layer(h_V_t, h_ESV_t, mask_V=mask_t))

      # Sampling step
      t_ = jnp.tile(t2, [1,1,h_V_stack[-1].shape[-1]])
      h_V_t = jnp.take_along_axis(h_V_stack[-1], t_, 1)[:,0]

      # sample
      logits = self.W_out(h_V_t) / temperature
      inputs = jnp.tile(jnp.arange(logits.shape[1])[None], [logits.shape[0], 1])
      key, sub_key = jax.random.split(key)
      S_t = jax.random.categorical(sub_key, logits)

      # updates
      t_ = jnp.tile(t2, [1,1,21])
      all_logits = scatter(all_logits, 1, t_, logits[:,None,:])
      h_S_ = self.W_s(S_t)
      t_ = jnp.tile(t2, [1,1,h_S_.shape[-1]])
      h_S = scatter(h_S, 1, t_, h_S_)
      S = scatter(S, 1, t[:,None], S_t)

    for i in range(N_nodes): step(i)

    output_dict = {"S": S, "logits": all_logits, "decoding_order": decoding_order}
    return output_dict