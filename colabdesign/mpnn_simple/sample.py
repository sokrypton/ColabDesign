import jax
import jax.numpy as jnp
import haiku as hk

import numpy as np
import itertools

from colabdesign.shared.prng import SafeKey
from .utils import gather_nodes, cat_neighbors_nodes, scatter, get_ar_mask

class mpnn_sample:
  def sample(self, X, mask, residue_idx, chain_idx, 
             decoding_order=None, temperature=0.1, bias=None):

    key = hk.next_rng_key()

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
    if decoding_order is None:
      key, sub_key = jax.random.split(key)
      decoding_order = jax.random.permutation(sub_key,jnp.arange(X.shape[1]))[None]
    ar_mask = get_ar_mask(decoding_order)

    mask_attend = jnp.take_along_axis(ar_mask, E_idx, 2)[...,None]
    mask_1D = mask.reshape([mask.shape[0], mask.shape[1], 1, 1])
    mask_bw = mask_1D * mask_attend
    mask_fw = mask_1D * (1. - mask_attend)
    
    h_EX_encoder = cat_neighbors_nodes(jnp.zeros_like(h_V), h_E, E_idx)
    h_EXV_encoder = cat_neighbors_nodes(h_V, h_EX_encoder, E_idx)
    h_EXV_encoder_fw = mask_fw * h_EXV_encoder

    def fn(x, t, key):
      
      # Hidden layers
      E_idx_t = E_idx[:,t][:,None]
      h_E_t = h_E[:,t][:,None]
      h_ES_t = cat_neighbors_nodes(x["h_S"], h_E_t, E_idx_t)
      h_EXV_encoder_t = h_EXV_encoder_fw[:,t][:,None] 
      mask_t = mask[:,t][:,None]
      mask_bw_t = mask_bw[:,t]

      def fn_layer(l, h_v):
        # Updated relational features for future states
        h_ESV_decoder_t = cat_neighbors_nodes(h_v[l], h_ES_t, E_idx_t)
        h_ESV_t = mask_bw_t * h_ESV_decoder_t + h_EXV_encoder_t

        # updates
        h_v_t = self.decoder_layers[l](h_v[l,:,t], h_ESV_t, mask_V=mask_t)[:,0]
        return h_v.at[l+1,:,t].set(h_v_t)

      for l in range(len(self.decoder_layers)):
        x["h_V_stack"] = fn_layer(l, x["h_V_stack"])

      # sample
      logits_t = self.W_out(x["h_V_stack"][-1,:,t])
      if bias is not None:
        logits_t = logits_t + bias[...,t,:]
      S_t_logits = logits_t/temperature + jax.random.gumbel(key, logits_t.shape)
      S_t = jax.nn.one_hot(S_t_logits.argmax(-1), 21)
      S_t_soft = jax.nn.softmax(S_t_logits,-1)
      S_t = jax.lax.stop_gradient(S_t - S_t_soft) + S_t_soft

      # updates
      x["h_S"] = x["h_S"].at[:,t].set(self.W_s(S_t))
      return x, (S_t, logits_t)
    
    # initial values
    N_batch, N_nodes = X.shape[0], X.shape[1]
    init = {"h_S":jnp.zeros_like(h_V),
            "h_V_stack":jnp.array([h_V] + [jnp.zeros_like(h_V)] * len(self.decoder_layers))}
    iters = {"t":decoding_order[0],
             "key":jax.random.split(key,N_nodes)}
    x, (S, logits) = jax.lax.scan(lambda x, xs: fn(x, **xs), init, iters)

    # put back in order
    i = jax.nn.one_hot(decoding_order[0], N_nodes).argmax(0)
    S,logits = S[i].swapaxes(0,1), logits[i].swapaxes(0,1)

    return {"S":S, "logits":logits, "decoding_order": decoding_order}