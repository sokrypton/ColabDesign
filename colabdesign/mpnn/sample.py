import jax
import jax.numpy as jnp
import haiku as hk
import numpy as np

from .utils import cat_neighbors_nodes, get_ar_mask

class mpnn_sample:
  def sample(self, I):

    key = hk.next_rng_key()

    # prepare node and edge embeddings
    E, E_idx = self.features(I)
    h_V = jnp.zeros((E.shape[0], E.shape[-1]))
    h_E = self.W_e(E)

    ##############
    # encoder
    ##############
    mask_attend = jnp.take_along_axis(I["mask"][:,None] * I["mask"][None,:], E_idx, 1)
    for layer in self.encoder_layers:
      h_V, h_E = layer(h_V, h_E, E_idx, I["mask"], mask_attend)

    if "decoding_order" not in I:
      key, sub_key = jax.random.split(key)
      randn = jax.random.uniform(sub_key, (I["X"].shape[0],))
      randn = jnp.where(I["mask"],randn,randn+1)
      if "fix_pos" in I: randn = randn.at[I["fix_pos"]].add(-1)
      I["decoding_order"] = randn.argsort()
    ar_mask = get_ar_mask(I["decoding_order"])

    mask_attend = jnp.take_along_axis(ar_mask, E_idx, 1)
    mask_1D = I["mask"][:,None]
    mask_bw = mask_1D * mask_attend
    mask_fw = mask_1D * (1 - mask_attend)
    
    h_EX_encoder = cat_neighbors_nodes(jnp.zeros_like(h_V), h_E, E_idx)
    h_EXV_encoder = cat_neighbors_nodes(h_V, h_EX_encoder, E_idx)
    h_EXV_encoder = mask_fw[...,None] * h_EXV_encoder

    def fn(x, t, key):      
      h_EXV_encoder_t = h_EXV_encoder[t,None] 
      E_idx_t         = E_idx[t,None]
      mask_t          = I["mask"][t,None]
      mask_bw_t       = mask_bw[t,None]      
      h_ES_t          = cat_neighbors_nodes(x["h_S"], h_E[t,None], E_idx_t)

      ##############
      # decoder
      ##############
      for l,layer in enumerate(self.decoder_layers):
        h_V = x["h_V"][l]
        h_ESV_decoder_t = cat_neighbors_nodes(h_V, h_ES_t, E_idx_t)
        h_ESV_t = mask_bw_t[...,None] * h_ESV_decoder_t + h_EXV_encoder_t
        h_V_t = layer(h_V[t], h_ESV_t, mask_V=mask_t)[0]
        # update
        x["h_V"] = x["h_V"].at[l+1,t].set(h_V_t)

      ##############
      # sample
      ##############
      logits_t = self.W_out(h_V_t)

      # add bias
      if "bias" in I:
        S_logits_t = logits_t + I["bias"][t]
      else:
        S_logits_t = logits_t

      # sample character
      S_logits_t = S_logits_t/I["temperature"] + jax.random.gumbel(key, S_logits_t.shape)
      S_t = jax.nn.one_hot(S_logits_t[:20].argmax(-1), 21)

      # backprop through sampling step
      S_soft_t = jnp.pad(jax.nn.softmax(S_logits_t[:20],-1),[0,1])
      S_t = jax.lax.stop_gradient(S_t - S_soft_t) + S_soft_t

      # update
      x["h_S"] = x["h_S"].at[t].set(self.W_s(S_t))
      return x, (S_t, logits_t)
    
    # initial values
    N_nodes = I["X"].shape[0]
    init = {"h_S":jnp.zeros_like(h_V),
            "h_V":jnp.array([h_V] + [jnp.zeros_like(h_V)] * len(self.decoder_layers))}
    iters = {"t":I["decoding_order"], "key":jax.random.split(key,N_nodes)}
    S, logits = jax.lax.scan(lambda x, xs: fn(x, **xs), init, iters)[1]

    # put back in order
    i = jax.nn.one_hot(I["decoding_order"], N_nodes).argmax(0)
    S,logits = S[i], logits[i]

    return {"S":S, "logits":logits, "decoding_order": I["decoding_order"]}