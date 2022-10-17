import jax
import jax.numpy as jnp
import numpy as np
import itertools

from .utils import gather_nodes, cat_neighbors_nodes, scatter, get_ar_mask

class mpnn_sample:
  def sample(self, key, X, randn, S_true,
         chain_mask, chain_idx, residue_idx,
         mask=None, temperature=1.0, omit_AAs_np=None,
         bias_AAs_np=None, chain_M_pos=None, omit_AA_mask=None,
         pssm_coef=None, pssm_bias=None, pssm_multi=None,
         pssm_log_odds_flag=None, pssm_log_odds_mask=None,
         pssm_bias_flag=None, bias_by_res=None):

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
    chain_mask = chain_mask * chain_M_pos * mask #update chain_M to include missing regions
    decoding_order = jnp.argsort((chain_mask+0.0001)*(jnp.abs(randn))) #[numbers will be smaller for places where chain_M = 0.0 and higher for places where chain_M = 1.0]

    ar_mask = get_ar_mask(decoding_order)

    mask_attend = jnp.take_along_axis(ar_mask, E_idx, 2)[...,None]
    mask_1D = mask.reshape([mask.shape[0], mask.shape[1], 1, 1])
    mask_bw = mask_1D * mask_attend
    mask_fw = mask_1D * (1. - mask_attend)

    N_batch, N_nodes = X.shape[0], X.shape[1]
    log_probs = jnp.zeros((N_batch, N_nodes, 21))
    all_probs = jnp.zeros((N_batch, N_nodes, 21))
    h_S = jnp.zeros_like(h_V,)
    S = jnp.zeros((N_batch, N_nodes), dtype=jnp.int32)
    h_V_stack = [h_V] + [jnp.zeros_like(h_V) for _ in range(len(self.decoder_layers))]
    constant = jnp.array(omit_AAs_np)
    constant_bias = jnp.array(bias_AAs_np)
    #chain_mask_combined = chain_mask*chain_M_pos 
    omit_AA_mask_flag = omit_AA_mask != None

    h_EX_encoder = cat_neighbors_nodes(jnp.zeros_like(h_S), h_E, E_idx)
    h_EXV_encoder = cat_neighbors_nodes(h_V, h_EX_encoder, E_idx)
    h_EXV_encoder_fw = mask_fw * h_EXV_encoder

    for t_ in range(N_nodes):
      t = decoding_order[:, t_]  # [B]
      chain_mask_gathered = jnp.take_along_axis(chain_mask, t[:,None], 1)  # [B]
      bias_by_res_gathered = jnp.take_along_axis(bias_by_res, jnp.tile(t[:,None, None], [1,1,21]), 1)[:,0,:]  # [B, 21]
      if jnp.equal(chain_mask_gathered, 0).all():
        S_t = jnp.take_along_axis(S_true, t[:,None], 1)
      else:
        # Hidden layers
        E_idx_t = jnp.take_along_axis(E_idx, jnp.tile(t[:,None,None], [1,1,E_idx.shape[-1]]), 1)
        h_E_t = jnp.take_along_axis(h_E, jnp.tile(t[:,None,None,None], [1,1,h_E.shape[-2], h_E.shape[-1]]), 1)
        h_ES_t = cat_neighbors_nodes(h_S, h_E_t, E_idx_t)
        h_EXV_encoder_t = jnp.take_along_axis(h_EXV_encoder_fw, jnp.tile(t[:,None,None,None], [1,1,h_EXV_encoder_fw.shape[-2], h_EXV_encoder_fw.shape[-1]]), 1)
        mask_t = jnp.take_along_axis(mask, t[:,None], 1)
        for l, layer in enumerate(self.decoder_layers):
          # Updated relational features for future states
          h_ESV_decoder_t = cat_neighbors_nodes(h_V_stack[l], h_ES_t, E_idx_t)
          h_V_t = jnp.take_along_axis(h_V_stack[l], jnp.tile(t[:,None,None], [1,1,h_V_stack[l].shape[-1]]), 1)
          h_ESV_t = jnp.take_along_axis(mask_bw, jnp.tile(t[:,None,None,None], [1,1,mask_bw.shape[-2], mask_bw.shape[-1]]), 1) * h_ESV_decoder_t + h_EXV_encoder_t
          h_V_stack[l+1] = scatter(h_V_stack[l+1], 1, jnp.tile(t[:,None,None], [1,1,h_V.shape[-1]]), layer(h_V_t, h_ESV_t, mask_V=mask_t))
        # Sampling step
        h_V_t = jnp.take_along_axis(h_V_stack[-1], jnp.tile(t[:,None,None], [1,1,h_V_stack[-1].shape[-1]]), 1)[:,0]
        logits = self.W_out(h_V_t) / temperature
        probs = jax.nn.softmax(logits-constant[None,:]*1e8+constant_bias[None,:]/temperature+bias_by_res_gathered/temperature, axis=-1)
        if pssm_bias_flag:
          pssm_coef_gathered = jnp.take_along_axis(pssm_coef, t[:,None], 1)[:,0]
          pssm_bias_gathered = jnp.take_along_axis(pssm_bias, jnp.tile(t[:,None,None], [1,1,pssm_bias.shape[-1]]), 1)[:,0]
          probs = (1-pssm_multi*pssm_coef_gathered[:,None])*probs + pssm_multi*pssm_coef_gathered[:,None]*pssm_bias_gathered
        if pssm_log_odds_flag:
          pssm_log_odds_mask_gathered = jnp.take_along_axis(pssm_log_odds_mask, jnp.tile(t[:,None, None], [1,1,pssm_log_odds_mask.shape[-1]]), 1)[:,0] #[B, 21]
          probs_masked = probs*pssm_log_odds_mask_gathered
          probs_masked += probs * 0.001
          probs = probs_masked/jnp.sum(probs_masked, axis=-1, keepdims=True) #[B, 21]
        if omit_AA_mask_flag:
          omit_AA_mask_gathered = jnp.take_along_axis(omit_AA_mask, jnp.tile(t[:,None, None], [1,1,omit_AA_mask.shape[-1]]), 1)[:,0] #[B, 21]
          probs_masked = probs*(1.0-omit_AA_mask_gathered)
          probs = probs_masked/jnp.sum(probs_masked, axis=-1, keepdims=True) #[B, 21]
        used_key = jax.random.split(key, probs.shape[0])
        input = jnp.tile(jnp.arange(probs.shape[1])[None], [probs.shape[0], 1])
        S_t = jax.vmap(lambda key, input, prob: jax.random.choice(key, input, p=prob),
                 in_axes=(0, 0, 0), out_axes=0)(used_key, input, probs)
        all_probs = scatter(all_probs, 1, jnp.tile(t[:,None,None], [1,1,21]),
          (chain_mask_gathered[:,:,None,]*probs[:,None,:]))
      S_true_gathered = jnp.take_along_axis(S_true, t[:,None], 1)
      S_t = (S_t*chain_mask_gathered+S_true_gathered*(1.0-chain_mask_gathered)).astype(int)
      temp1 = self.W_s(S_t)
      h_S = scatter(h_S, 1, jnp.tile(t[:,None,None], [1,1,temp1.shape[-1]]), temp1)
      S = scatter(S, 1, t[:,None], S_t)
    output_dict = {"S": S, "probs": all_probs, "decoding_order": decoding_order}
    return output_dict


  def tied_sample(self, key, X, randn, S_true,
          chain_mask, chain_idx, residue_idx,
          mask=None, temperature=1.0, omit_AAs_np=None,
          bias_AAs_np=None, chain_M_pos=None, omit_AA_mask=None,
          pssm_coef=None, pssm_bias=None, pssm_multi=None,
          pssm_log_odds_flag=None, pssm_log_odds_mask=None,
          pssm_bias_flag=None, tied_pos=None, tied_beta=None,
          bias_by_res=None):

    # Prepare node and edge embeddings
    E, E_idx = self.features(X, mask, residue_idx, chain_idx)
    h_V = jnp.zeros((E.shape[0], E.shape[1], E.shape[-1]))
    h_E = self.W_e(E)
    # Encoder is unmasked self-attention 
    mask_attend = gather_nodes(mask[...,None],E_idx)[...,0]
    mask_attend = mask[...,None] * mask_attend
    for layer in self.encoder_layers:
      h_V, h_E = layer(h_V, h_E, E_idx, mask, mask_attend)

    # Decoder uses masked self-attention
    chain_mask = chain_mask*chain_M_pos*mask #update chain_M to include missing regions
    decoding_order = jnp.argsort((chain_mask+0.0001)*(jnp.abs(randn))) #[numbers will be smaller for places where chain_M = 0.0 and higher for places where chain_M = 1.0]

    new_decoding_order = []
    for t_dec in list(decoding_order[0]):
      if t_dec not in list(itertools.chain(*new_decoding_order)):
        list_a = [item for item in tied_pos if t_dec in item]
        if list_a:
          new_decoding_order.append(list_a[0])
        else:
          new_decoding_order.append([t_dec])
    decoding_order = jnp.tile(jnp.array(list(itertools.chain(*new_decoding_order)))[None,], [X.shape[0], 1])

    ar_mask = get_ar_mask(decoding_order)

    mask_attend = jnp.take_along_axis(ar_mask, E_idx, 2)[...,None]
    mask_1D = mask.reshape([mask.shape[0], mask.shape[1], 1, 1])
    mask_bw = mask_1D * mask_attend
    mask_fw = mask_1D * (1. - mask_attend)

    N_batch, N_nodes = X.shape[0], X.shape[1]
    log_probs = jnp.zeros((N_batch, N_nodes, 21))
    all_probs = jnp.zeros((N_batch, N_nodes, 21))
    h_S = jnp.zeros_like(h_V)
    S = jnp.zeros((N_batch, N_nodes), dtype=jnp.int32)
    h_V_stack = [h_V] + [jnp.zeros_like(h_V) for _ in range(len(self.decoder_layers))]
    constant = jnp.array(omit_AAs_np)
    constant_bias = jnp.array(bias_AAs_np)
    omit_AA_mask_flag = omit_AA_mask != None

    h_EX_encoder = cat_neighbors_nodes(jnp.zeros_like(h_S), h_E, E_idx)
    h_EXV_encoder = cat_neighbors_nodes(h_V, h_EX_encoder, E_idx)
    h_EXV_encoder_fw = mask_fw * h_EXV_encoder
    for t_list in new_decoding_order:
      logits = 0.0
      logit_list = []
      done_flag = False
      for t in t_list:
        if (chain_mask[:,t]==0).all():
          S_t = S_true[:,t]
          for t in t_list:
            h_S[:,t,:] = self.W_s(S_t)
            S[:,t] = S_t
          done_flag = True
          break
        else:
          E_idx_t = E_idx[:,t:t+1,:]
          h_E_t = h_E[:,t:t+1,:,:]
          h_ES_t = cat_neighbors_nodes(h_S, h_E_t, E_idx_t)
          h_EXV_encoder_t = h_EXV_encoder_fw[:,t:t+1,:,:]
          mask_t = mask[:,t:t+1]
          for l, layer in enumerate(self.decoder_layers):
            h_ESV_decoder_t = cat_neighbors_nodes(h_V_stack[l], h_ES_t, E_idx_t)
            h_V_t = h_V_stack[l][:,t:t+1,:]
            h_ESV_t = mask_bw[:,t:t+1,:,:] * h_ESV_decoder_t + h_EXV_encoder_t
            h_V_stack[l+1] = h_V_stack[l+1].at[:,t,:].set(layer(h_V_t, h_ESV_t, mask_V=mask_t).squeeze(1))
            # h_V_stack[l+1][:,t,:] = layer(h_V_t, h_ESV_t, mask_V=mask_t).squeeze(1)
          h_V_t = h_V_stack[-1][:,t,:]
          logit_list.append((self.W_out(h_V_t) / temperature)/len(t_list))
          logits += tied_beta[t]*(self.W_out(h_V_t) / temperature)/len(t_list)
      if done_flag:
        pass
      else:
        bias_by_res_gathered = bias_by_res[:,t,:] #[B, 21]
        probs = jax.nn.softmax(logits-constant[None,:]*1e8+constant_bias[None,:]/temperature+bias_by_res_gathered/temperature, axis=-1)
        if pssm_bias_flag:
          pssm_coef_gathered = pssm_coef[:,t]
          pssm_bias_gathered = pssm_bias[:,t]
          probs = (1-pssm_multi*pssm_coef_gathered[:,None])*probs + pssm_multi*pssm_coef_gathered[:,None]*pssm_bias_gathered
        if pssm_log_odds_flag:
          pssm_log_odds_mask_gathered = pssm_log_odds_mask[:,t]
          probs_masked = probs*pssm_log_odds_mask_gathered
          probs_masked += probs * 0.001
          probs = probs_masked/jnp.sum(probs_masked, aixs=-1, keepdims=True) #[B, 21]
        if omit_AA_mask_flag:
          omit_AA_mask_gathered = omit_AA_mask[:,t]
          probs_masked = probs*(1.0-omit_AA_mask_gathered)
          probs = probs_masked/jnp.sum(probs_masked, axis=-1, keepdims=True) #[B, 21]
        
        used_key = jax.random.split(key, probs.shape[0])
        input = jnp.tile(jnp.arange(probs.shape[1])[None], [probs.shape[0], 1])
        S_t_repeat = jax.vmap(lambda key, input, prob: jax.random.choice(key, input, p=prob),
                    in_axes=(0, 0, 0), out_axes=0)(used_key, input, probs)

        for t in t_list:
          h_S = h_S.at[:,t,:].set(self.W_s(S_t_repeat))
          S = S.at[:,t].set(S_t_repeat)
          all_probs = all_probs.at[:,t,:].set(probs)
    output_dict = {"S": S, "probs": all_probs, "decoding_order": decoding_order}
    return output_dict