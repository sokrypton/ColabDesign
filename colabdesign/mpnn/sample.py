import jax
import jax.numpy as jnp
import haiku as hk
import numpy as np

from .utils import cat_neighbors_nodes, get_ar_mask

class mpnn_sample:
  def sample(self, I, tied_lengths=False):

    key = hk.next_rng_key()
    N_nodes = I["X"].shape[0]

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

    ######################
    # get decoding order
    ######################
    key, sub_key = jax.random.split(key)

    # get decoding order
    if "decoding_order" in I:
      decoding_order = I["decoding_order"]
      if decoding_order.ndim == 1:
        decoding_order = decoding_order[:,None]

    else:
      randn = jax.random.uniform(key, (N_nodes,))    
      randn = jnp.where(I["mask"],randn,randn + 1)
      if "fix_pos" in I: randn = randn.at[I["fix_pos"]].add(-1)
      if tied_lengths:
        lengths = I.get("lengths",jnp.array([N_nodes]))    
        copies = lengths.shape[0]
        decoding_order_tied = randn.reshape(copies,-1).mean(0).argsort()
        decoding_order = jnp.arange(N_nodes).reshape(copies,-1).T[decoding_order_tied]
      else:
        decoding_order = randn.argsort(-1)[:,None]
    
    # get autoregressive mask
    if "ar_mask" in I:
      ar_mask = I["ar_mask"]
    else:
      ar_mask = get_ar_mask(decoding_order)
    
    mask_attend = jnp.take_along_axis(ar_mask, E_idx, 1)
    mask_1D = I["mask"][:,None]
    mask_bw = mask_1D * mask_attend
    mask_fw = mask_1D * (1 - mask_attend)
    
    h_EX_encoder = cat_neighbors_nodes(jnp.zeros_like(h_V), h_E, E_idx)
    h_EXV_encoder = cat_neighbors_nodes(h_V, h_EX_encoder, E_idx)
    h_EXV_encoder = mask_fw[...,None] * h_EXV_encoder

    def fn(x, t, key, tied=True, sample=True):
      h_EXV_encoder_t = h_EXV_encoder[t] 
      E_idx_t         = E_idx[t]
      mask_t          = I["mask"][t]
      mask_bw_t       = mask_bw[t]      
      h_ES_t          = cat_neighbors_nodes(x["h_S"], h_E[t], E_idx_t)

      ##############
      # decoder
      ##############
      for l,layer in enumerate(self.decoder_layers):
        h_V = x["h_V"][l]
        h_ESV_decoder_t = cat_neighbors_nodes(h_V, h_ES_t, E_idx_t)
        h_ESV_t = mask_bw_t[...,None] * h_ESV_decoder_t + h_EXV_encoder_t
        h_V_t = layer(h_V[t], h_ESV_t, mask_V=mask_t)
        # update
        x["h_V"] = x["h_V"].at[l+1,t].set(h_V_t)

      logits_t = self.W_out(h_V_t)
      if tied: logits_t = logits_t.mean(0,keepdims=True)
      x["logits"] = x["logits"].at[t].set(logits_t)

      ##############
      # sample
      ##############
      if sample:
        # add bias
        if "bias" in I:
          S_logits_t = logits_t + I["bias"][t]
        else:
          S_logits_t = logits_t

        if tied: S_logits_t = S_logits_t.mean(0,keepdims=True)

        # sample character
        S_logits_t = S_logits_t/I["temperature"] + jax.random.gumbel(key, S_logits_t.shape)
        S_t = jax.nn.one_hot(S_logits_t[...,:20].argmax(-1), 21)

        # backprop through sampling step
        S_soft_t = jnp.pad(jax.nn.softmax(S_logits_t[...,:20],-1),[[0,0],[0,1]])
        S_t = jax.lax.stop_gradient(S_t - S_soft_t) + S_soft_t

        # update
        x["h_S"] = x["h_S"].at[t].set(self.W_s(S_t))
        x["S"]   = x["S"].at[t].set(S_t)
      return x
    
    # initial values
    O = {"h_S":    jnp.zeros_like(h_V),
         "h_V":    jnp.array([h_V] + [jnp.zeros_like(h_V)] * len(self.decoder_layers)),
         "S":      jnp.zeros((N_nodes,21)),
         "logits": jnp.zeros((N_nodes,21))}

    # scan over decoding order
    iters = {"t":decoding_order,
             "key":jax.random.split(key,decoding_order.shape[0])}
    O = hk.scan(lambda x, xs: (fn(x,**xs),None), O, iters)[0]        

    return {"S":O["S"], "logits":O["logits"], "decoding_order":decoding_order}
