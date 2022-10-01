import functools
import haiku as hk
import jax
import jax.numpy as jnp
import joblib
import numpy as np
import itertools

from colabdesign.mpnn.lib.prng import SafeKey

Gelu = functools.partial(jax.nn.gelu, approximate=False)

# The following gather functions
def gather_edges(edges, neighbor_idx):
    # Features [B,N,N,C] at Neighbor indices [B,N,K] => Neighbor features [B,N,K,C]
    neighbors = jnp.tile(jnp.expand_dims(neighbor_idx, -1), [1, 1, 1, edges.shape[-1]])
    edge_features = jnp.take_along_axis(edges, neighbors, 2)
    return edge_features

def gather_nodes(nodes, neighbor_idx):
    # Features [B,N,C] at Neighbor indices [B,N,K] => [B,N,K,C]
    # Flatten and expand indices per batch [B,N,K] => [B,NK] => [B,NK,C]
    neighbors_flat = neighbor_idx.reshape([neighbor_idx.shape[0], -1])
    neighbors_flat = jnp.tile(jnp.expand_dims(neighbors_flat, -1),[1, 1, nodes.shape[2]])
    # Gather and re-pack
    neighbor_features = jnp.take_along_axis(nodes, neighbors_flat, 1)
    neighbor_features = neighbor_features.reshape(list(neighbor_idx.shape[:3]) + [-1])
    return neighbor_features

def gather_nodes_t(nodes, neighbor_idx):
    # Features [B,N,C] at Neighbor index [B,K] => Neighbor features[B,K,C]
    idx_flat = jnp.tile(jnp.expand_dims(neighbor_idx, -1),[1, 1, nodes.shape[2]])
    neighbor_features = jnp.take_along_axis(nodes, idx_flat, 1)
    return neighbor_features

def cat_neighbors_nodes(h_nodes, h_neighbors, E_idx):
    h_nodes = gather_nodes(h_nodes, E_idx)
    h_nn = jnp.concatenate([h_neighbors, h_nodes], -1)
    return h_nn

def scatter(input, dim, index, src):
    idx = jnp.indices(index.shape)
    dim_num = idx.shape[0]  # dimension of the dim
    a = idx.at[dim].set(index)
    idx = jnp.moveaxis(idx, 0, -1).reshape((-1, dim_num))
    a = jnp.moveaxis(a, 0, -1).reshape((-1, dim_num))
    return input.at[tuple(a.T)].set(src[tuple(idx.T)])


class dropout_cust(hk.Module):
    def __init__(self, rate) -> None:
        super().__init__()
        self.rate = rate
        self.safe_key = SafeKey(hk.next_rng_key())
    
    def __call__(self, x):
        self.safe_key, use_key = self.safe_key.split()
        return hk.dropout(use_key.get(), self.rate, x)


class EncLayer(hk.Module):
    def __init__(self, num_hidden,
                 num_in, dropout=0.1,
                 num_heads=None, scale=30,
                 name=None):
        super(EncLayer, self).__init__()
        self.num_hidden = num_hidden
        self.num_in = num_in
        self.scale = scale

        self.safe_key = SafeKey(hk.next_rng_key())

        self.dropout1 = dropout_cust(dropout)
        self.dropout2 = dropout_cust(dropout)
        self.dropout3 = dropout_cust(dropout)
        self.norm1 = hk.LayerNorm(-1, create_scale=True, create_offset=True,
                                  name=name + '_norm1')
        self.norm2 = hk.LayerNorm(-1, create_scale=True, create_offset=True,
                                  name=name + '_norm2')
        self.norm3 = hk.LayerNorm(-1, create_scale=True, create_offset=True,
                                  name=name + '_norm3')

        self.W1 = hk.Linear(num_hidden, with_bias=True, name=name + '_W1')
        self.W2 = hk.Linear(num_hidden, with_bias=True, name=name + '_W2')
        self.W3 = hk.Linear(num_hidden, with_bias=True, name=name + '_W3')
        self.W11 = hk.Linear(num_hidden, with_bias=True, name=name + '_W11')
        self.W12 = hk.Linear(num_hidden, with_bias=True, name=name + '_W12')
        self.W13 = hk.Linear(num_hidden, with_bias=True, name=name + '_W13')
        self.act = Gelu
        self.dense = PositionWiseFeedForward(num_hidden, num_hidden * 4,
                                             name=name + '_dense')

    def __call__(self, h_V, h_E, E_idx,
                mask_V=None, mask_attend=None):
        """ Parallel computation of full transformer layer """

        h_EV = cat_neighbors_nodes(h_V, h_E, E_idx)
        h_V_expand = jnp.tile(jnp.expand_dims(h_V, -2),[1, 1, h_EV.shape[-2], 1])
        h_EV = jnp.concatenate([h_V_expand, h_EV], -1)

        h_message = self.W3(self.act(self.W2(self.act(self.W1(h_EV)))))
        if mask_attend is not None:
            h_message = jnp.expand_dims(mask_attend, -1)* h_message
        dh = jnp.sum(h_message, -2) / self.scale
        h_V = self.norm1(h_V + self.dropout1(dh))

        dh = self.dense(h_V)
        h_V = self.norm2(h_V + self.dropout2(dh))
        if mask_V is not None:
            mask_V = jnp.expand_dims(mask_V, -1)
            h_V = mask_V * h_V

        h_EV = cat_neighbors_nodes(h_V, h_E, E_idx)
        h_V_expand = jnp.tile(jnp.expand_dims(h_V, -2),[1, 1, h_EV.shape[-2], 1])
        h_EV = jnp.concatenate([h_V_expand, h_EV], -1)
        h_message = self.W13(self.act(self.W12(self.act(self.W11(h_EV)))))
        h_E = self.norm3(h_E + self.dropout3(h_message))
        return h_V, h_E


class DecLayer(hk.Module):
    def __init__(self, num_hidden, num_in,
                 dropout=0.1, num_heads=None,
                 scale=30, name=None):
        super(DecLayer, self).__init__()
        self.num_hidden = num_hidden
        self.num_in = num_in
        self.scale = scale
        self.dropout1 = dropout_cust(dropout)
        self.dropout2 = dropout_cust(dropout)
        self.norm1 = hk.LayerNorm(-1, create_scale=True, create_offset=True,
                                  name=name + '_norm1')
        self.norm2 = hk.LayerNorm(-1, create_scale=True, create_offset=True,
                                  name=name + '_norm2')

        self.W1 = hk.Linear(num_hidden, with_bias=True, name=name + '_W1')
        self.W2 = hk.Linear(num_hidden, with_bias=True, name=name + '_W2')
        self.W3 = hk.Linear(num_hidden, with_bias=True, name=name + '_W3')
        self.act = Gelu
        self.dense = PositionWiseFeedForward(num_hidden, num_hidden * 4,
                                             name=name + '_dense')

    def __call__(self, h_V, h_E,
                 mask_V=None, mask_attend=None):
        """ Parallel computation of full transformer layer """

        # Concatenate h_V_i to h_E_ij
        h_V_expand = jnp.tile(jnp.expand_dims(h_V, -2),[1, 1, h_E.shape[-2], 1])
        h_EV = jnp.concatenate([h_V_expand, h_E], -1)

        h_message = self.W3(self.act(self.W2(self.act(self.W1(h_EV)))))
        if mask_attend is not None:
            h_message = jnp.expand_dims(mask_attend, -1) * h_message
        dh = jnp.sum(h_message, -2) / self.scale

        h_V = self.norm1(h_V + self.dropout1(dh))

        # Position-wise feedforward
        dh = self.dense(h_V)
        h_V = self.norm2(h_V + self.dropout2(dh))

        if mask_V is not None:
            mask_V = jnp.expand_dims(mask_V, -1)
            h_V = mask_V * h_V
        return h_V 


class PositionWiseFeedForward(hk.Module):
    def __init__(self, num_hidden, num_ff, name=None):
        super(PositionWiseFeedForward, self).__init__()
        self.W_in = hk.Linear(num_ff, with_bias=True, name=name + '_W_in')
        self.W_out = hk.Linear(num_hidden, with_bias=True, name=name + '_W_out')
        self.act = Gelu
    def __call__(self, h_V):
        h = self.act(self.W_in(h_V), approximate=False)
        h = self.W_out(h)
        return h


class PositionalEncodings(hk.Module):
    def __init__(self, num_embeddings, max_relative_feature=32):
        super(PositionalEncodings, self).__init__()
        self.num_embeddings = num_embeddings
        self.max_relative_feature = max_relative_feature
        self.linear = hk.Linear(num_embeddings, name='embedding_linear')

    def __call__(self, offset, mask):
        d = jnp.clip(offset + self.max_relative_feature, 0, 2*self.max_relative_feature) * mask + \
            (1 - mask) * (2*self.max_relative_feature + 1)
        d_onehot = jax.nn.one_hot(d, 2*self.max_relative_feature + 1 + 1)
        E = self.linear(jax.lax.convert_element_type(d_onehot, jnp.float32))
        return E


class RunModel:
    def __init__(self, config) -> None:
        self.config = config

        def _forward_score(input):
            model = ProteinMPNN(**self.config)
            return model(**input)
        self.score = jax.jit(hk.transform(_forward_score).apply)
        self.init_score = jax.jit(hk.transform(_forward_score).init)

        def _forward_sample(input):
            model = ProteinMPNN(**self.config)
            return model.sample(**input)
        self.sample = hk.transform(_forward_sample).apply
        self.init_sample = hk.transform(_forward_sample).init

        def _forward_tsample(input):
            model = ProteinMPNN(**self.config)
            return model.tied_sample(**input)
        self.tied_sample = hk.transform(_forward_tsample).apply
        self.init_tsample = hk.transform(_forward_tsample).init

    def load_params(self, path):
        self.params = joblib.load(path)


class ProteinFeatures(hk.Module):
    def __init__(self, edge_features, node_features,
                 num_positional_embeddings=16,
                 num_rbf=16, top_k=30,
                 augment_eps=0., num_chain_embeddings=16):
        """ Extract protein features """
        super(ProteinFeatures, self).__init__()
        self.edge_features = edge_features
        self.node_features = node_features
        self.top_k = top_k
        self.augment_eps = augment_eps 
        self.num_rbf = num_rbf
        self.num_positional_embeddings = num_positional_embeddings

        self.embeddings = PositionalEncodings(num_positional_embeddings)
        node_in, edge_in = 6, num_positional_embeddings + num_rbf*25
        self.edge_embedding = hk.Linear(edge_features, with_bias=False, name='edge_embedding')
        self.norm_edges = hk.LayerNorm(-1, create_scale=True, create_offset=True, name='norm_edges')

        self.safe_key = SafeKey(hk.next_rng_key())

    def _dist(self, X, mask, eps=1E-6):
        mask_2D = jnp.expand_dims(mask, 1) * jnp.expand_dims(mask, 2)
        dX = jnp.expand_dims(X, 1) - jnp.expand_dims(X, 2)
        D = mask_2D * jnp.sqrt(jnp.sum(dX**2, 3) + eps)
        D_max = jnp.max(D, -1, keepdims=True)
        D_adjust = D + (1. - mask_2D) * D_max
        sampled_top_k = self.top_k
        D_neighbors, E_idx = jax.lax.approx_min_k(D_adjust,
                                                  np.minimum(self.top_k, X.shape[1]),
                                                  reduction_dimension=-1)
        return D_neighbors, E_idx

    def _rbf(self, D):
        D_min, D_max, D_count = 2., 22., self.num_rbf
        D_mu = jnp.linspace(D_min, D_max, D_count)
        D_mu = D_mu.reshape([1,1,1,-1])
        D_sigma = (D_max - D_min) / D_count
        D_expand = jnp.expand_dims(D, -1)
        RBF = jnp.exp(-((D_expand - D_mu) / D_sigma)**2)
        return RBF

    def _get_rbf(self, A, B, E_idx):
        D_A_B = jnp.sqrt(jnp.sum((A[:,:,None,:] - B[:,None,:,:])**2, -1) + 1e-6) #[B, L, L]
        D_A_B_neighbors = gather_edges(D_A_B[:,:,:,None], E_idx)[:,:,:,0] #[B,L,K]
        RBF_A_B = self._rbf(D_A_B_neighbors)
        return RBF_A_B

    def __call__(self, X, mask, residue_idx, chain_labels):
        if self.augment_eps > 0:
            self.safe_key, use_key = self.safe_key.split()
            X = X + self.augment_eps * jax.random.normal(use_key, X.shape)
        
        b = X[:,:,1,:] - X[:,:,0,:]
        c = X[:,:,2,:] - X[:,:,1,:]
        a = jnp.cross(b, c)
        Cb = -0.58273431*a + 0.56802827*b - 0.54067466*c + X[:,:,1,:]
        Ca = X[:,:,1,:]
        N = X[:,:,0,:]
        C = X[:,:,2,:]
        O = X[:,:,3,:]

        D_neighbors, E_idx = self._dist(Ca, mask)

        RBF_all = []
        RBF_all.append(self._rbf(D_neighbors)) #Ca-Ca
        RBF_all.append(self._get_rbf(N, N, E_idx)) #N-N
        RBF_all.append(self._get_rbf(C, C, E_idx)) #C-C
        RBF_all.append(self._get_rbf(O, O, E_idx)) #O-O
        RBF_all.append(self._get_rbf(Cb, Cb, E_idx)) #Cb-Cb
        RBF_all.append(self._get_rbf(Ca, N, E_idx)) #Ca-N
        RBF_all.append(self._get_rbf(Ca, C, E_idx)) #Ca-C
        RBF_all.append(self._get_rbf(Ca, O, E_idx)) #Ca-O
        RBF_all.append(self._get_rbf(Ca, Cb, E_idx)) #Ca-Cb
        RBF_all.append(self._get_rbf(N, C, E_idx)) #N-C
        RBF_all.append(self._get_rbf(N, O, E_idx)) #N-O
        RBF_all.append(self._get_rbf(N, Cb, E_idx)) #N-Cb
        RBF_all.append(self._get_rbf(Cb, C, E_idx)) #Cb-C
        RBF_all.append(self._get_rbf(Cb, O, E_idx)) #Cb-O
        RBF_all.append(self._get_rbf(O, C, E_idx)) #O-C
        RBF_all.append(self._get_rbf(N, Ca, E_idx)) #N-Ca
        RBF_all.append(self._get_rbf(C, Ca, E_idx)) #C-Ca
        RBF_all.append(self._get_rbf(O, Ca, E_idx)) #O-Ca
        RBF_all.append(self._get_rbf(Cb, Ca, E_idx)) #Cb-Ca
        RBF_all.append(self._get_rbf(C, N, E_idx)) #C-N
        RBF_all.append(self._get_rbf(O, N, E_idx)) #O-N
        RBF_all.append(self._get_rbf(Cb, N, E_idx)) #Cb-N
        RBF_all.append(self._get_rbf(C, Cb, E_idx)) #C-Cb
        RBF_all.append(self._get_rbf(O, Cb, E_idx)) #O-Cb
        RBF_all.append(self._get_rbf(C, O, E_idx)) #C-O
        RBF_all = jnp.concatenate(tuple(RBF_all), axis=-1)

        offset = residue_idx[:,:,None] - residue_idx[:,None,:]
        offset = gather_edges(offset[:,:,:,None], E_idx)[:,:,:,0] #[B, L, K]

        d_chains = (chain_labels[:, :, None] - chain_labels[:,None,:])==0
        d_chains = jax.lax.convert_element_type(d_chains, jnp.int64)
        E_chains = gather_edges(d_chains[:,:,:,None], E_idx)[:,:,:,0]
        E_positional = self.embeddings(jax.lax.convert_element_type(offset, jnp.int64), E_chains)
        E = jnp.concatenate((E_positional, RBF_all), -1)
        E = self.edge_embedding(E)
        E = self.norm_edges(E)
        return E, E_idx 


class EmbedToken(hk.Module):
    def __init__(self, vocab_size, embed_dim):
        super().__init__()
        self.vocab_size = vocab_size
        self.embed_dim = embed_dim
        self.w_init = hk.initializers.TruncatedNormal()

    @property
    def embeddings(self):
        return hk.get_parameter("W_s",
                                [self.vocab_size, self.embed_dim],
                                init=self.w_init)

    def __call__(self, arr):
        one_hot = jax.nn.one_hot(arr, self.vocab_size)
        return jnp.tensordot(one_hot, self.embeddings, 1)


class ProteinMPNN(hk.Module):
    def __init__(self, num_letters, node_features, 
                 edge_features, hidden_dim,
                 num_encoder_layers=3, num_decoder_layers=3,
                 vocab=21, k_neighbors=64,
                 augment_eps=0.05, dropout=0.1):
        super(ProteinMPNN, self).__init__()

        # Hyperparameters
        self.node_features = node_features
        self.edge_features = edge_features
        self.hidden_dim = hidden_dim

        # Featurization layers
        self.features = ProteinFeatures(node_features,
                                        edge_features,
                                        top_k=k_neighbors,
                                        augment_eps=augment_eps)

        self.W_e = hk.Linear(hidden_dim, with_bias=True, name='W_e')
        self.W_s = EmbedToken(vocab_size=vocab,
                              embed_dim=hidden_dim)

        # Encoder layers
        self.encoder_layers = [
            EncLayer(hidden_dim, hidden_dim*2, dropout=dropout, name='enc' + str(i))
            for i in range(num_encoder_layers)
        ]

        # Decoder layers
        self.decoder_layers = [
            DecLayer(hidden_dim, hidden_dim*3, dropout=dropout, name='dec' + str(i))
            for i in range(num_decoder_layers)
        ]
        self.W_out = hk.Linear(num_letters, with_bias=True, name='W_out')

    def __call__(self, X, S,
                 mask, chain_M,
                 residue_idx, chain_encoding_all,
                 randn, use_input_decoding_order=False,
                 decoding_order=None):
        """ Graph-conditioned sequence model """
        # Prepare node and edge embeddings
        E, E_idx = self.features(X, mask, residue_idx, chain_encoding_all)
        h_V = jnp.zeros((E.shape[0], E.shape[1], E.shape[-1]))
        h_E = self.W_e(E)

        # Encoder is unmasked self-attention
        mask_attend = gather_nodes(jnp.expand_dims(mask, -1),  E_idx).squeeze(-1)
        mask_attend = jnp.expand_dims(mask, -1) * mask_attend
        for layer in self.encoder_layers:
            h_V, h_E = layer(h_V, h_E, E_idx, mask, mask_attend)

        # Concatenate sequence embeddings for autoregressive decoder
        h_S = self.W_s(S)
        h_ES = cat_neighbors_nodes(h_S, h_E, E_idx)

        # Build encoder embeddings
        h_EX_encoder = cat_neighbors_nodes(jnp.zeros_like(h_S), h_E, E_idx)
        h_EXV_encoder = cat_neighbors_nodes(h_V, h_EX_encoder, E_idx)

        # update chain_M to include missing regions
        chain_M = chain_M * mask
        if not use_input_decoding_order:
            #[numbers will be smaller for places where chain_M = 0.0 and higher for places where chain_M = 1.0]
            decoding_order = jnp.argsort((chain_M+0.0001)*(jnp.abs(randn)))
        mask_size = E_idx.shape[1]
        permutation_matrix_reverse = jax.lax.convert_element_type(jax.nn.one_hot(decoding_order, mask_size),
                                                                  jnp.float32)
        order_mask_backward = jnp.einsum('ij, biq, bjp->bqp',
                                         (1 - jnp.triu(jnp.ones([mask_size, mask_size]))),
                                         permutation_matrix_reverse,
                                         permutation_matrix_reverse)
                      
        mask_attend = jnp.expand_dims(jnp.take_along_axis(order_mask_backward, E_idx, 2), -1)
        mask_1D = mask.reshape([mask.shape[0], mask.shape[1], 1, 1])
        mask_bw = mask_1D * mask_attend
        mask_fw = mask_1D * (1. - mask_attend)

        h_EXV_encoder_fw = mask_fw * h_EXV_encoder
        for layer in self.decoder_layers:
            # Masked positions attend to encoder information, unmasked see. 
            h_ESV = cat_neighbors_nodes(h_V, h_ES, E_idx)
            h_ESV = mask_bw * h_ESV + h_EXV_encoder_fw
            h_V = layer(h_V, h_ESV, mask)

        logits = self.W_out(h_V)
        log_probs = jax.nn.log_softmax(logits, axis=-1)
        return logits, log_probs
    
    def sample(self, key, X, randn, S_true,
               chain_mask, chain_encoding_all, residue_idx,
               mask=None, temperature=1.0, omit_AAs_np=None,
               bias_AAs_np=None, chain_M_pos=None, omit_AA_mask=None,
               pssm_coef=None, pssm_bias=None, pssm_multi=None,
               pssm_log_odds_flag=None, pssm_log_odds_mask=None,
               pssm_bias_flag=None, bias_by_res=None):

        # Prepare node and edge embeddings
        E, E_idx = self.features(X, mask, residue_idx, chain_encoding_all)
        h_V = jnp.zeros((E.shape[0], E.shape[1], E.shape[-1]))
        h_E = self.W_e(E)

        # Encoder is unmasked self-attention
        mask_attend = gather_nodes(jnp.expand_dims(mask, -1),  E_idx).squeeze(-1)
        mask_attend = jnp.expand_dims(mask, -1) * mask_attend
        for layer in self.encoder_layers:
            h_V, h_E = layer(h_V, h_E, E_idx, mask, mask_attend)

        # Decoder uses masked self-attention
        chain_mask = chain_mask*chain_M_pos*mask #update chain_M to include missing regions
        decoding_order = jnp.argsort((chain_mask+0.0001)*(jnp.abs(randn))) #[numbers will be smaller for places where chain_M = 0.0 and higher for places where chain_M = 1.0]
        mask_size = E_idx.shape[1]
        permutation_matrix_reverse = jax.lax.convert_element_type(jax.nn.one_hot(decoding_order, mask_size),
                                                                  jnp.float32)
        order_mask_backward = jnp.einsum('ij, biq, bjp->bqp',
                                         (1 - jnp.triu(jnp.ones([mask_size, mask_size]))),
                                         permutation_matrix_reverse,
                                         permutation_matrix_reverse)

        mask_attend = jnp.expand_dims(jnp.take_along_axis(order_mask_backward, E_idx, 2), -1)
        mask_1D = mask.reshape([mask.shape[0], mask.shape[1], 1, 1])
        mask_bw = mask_1D * mask_attend
        mask_fw = mask_1D * (1. - mask_attend)

        N_batch, N_nodes = X.shape[0], X.shape[1]
        log_probs = jnp.zeros((N_batch, N_nodes, 21))
        all_probs = jnp.zeros((N_batch, N_nodes, 21), dtype=jnp.float32)
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
            bias_by_res_gathered = jnp.take_along_axis(bias_by_res, t[:,None, None].tile([1,1,21]), 1)[:,0,:]  # [B, 21]
            if jnp.equal(chain_mask_gathered, 0).all():
                S_t = jnp.take_along_axis(S_true, t[:,None], 1)
            else:
                # Hidden layers
                E_idx_t = jnp.take_along_axis(E_idx, t[:,None,None].tile([1,1,E_idx.shape[-1]]), 1)
                h_E_t = jnp.take_along_axis(h_E, t[:,None,None,None].tile([1,1,h_E.shape[-2], h_E.shape[-1]]), 1)
                h_ES_t = cat_neighbors_nodes(h_S, h_E_t, E_idx_t)
                h_EXV_encoder_t = jnp.take_along_axis(h_EXV_encoder_fw, t[:,None,None,None].tile([1,1,h_EXV_encoder_fw.shape[-2], h_EXV_encoder_fw.shape[-1]]), 1)
                mask_t = jnp.take_along_axis(mask, t[:,None], 1)
                for l, layer in enumerate(self.decoder_layers):
                    # Updated relational features for future states
                    h_ESV_decoder_t = cat_neighbors_nodes(h_V_stack[l], h_ES_t, E_idx_t)
                    h_V_t = jnp.take_along_axis(h_V_stack[l], t[:,None,None].tile([1,1,h_V_stack[l].shape[-1]]), 1)
                    h_ESV_t = jnp.take_along_axis(mask_bw, t[:,None,None,None].tile([1,1,mask_bw.shape[-2], mask_bw.shape[-1]]), 1) * h_ESV_decoder_t + h_EXV_encoder_t
                    h_V_stack[l+1] = scatter(h_V_stack[l+1], 1, t[:,None,None].tile([1,1,h_V.shape[-1]]), layer(h_V_t, h_ESV_t, mask_V=mask_t))
                # Sampling step
                h_V_t = jnp.take_along_axis(h_V_stack[-1], t[:,None,None].tile([1,1,h_V_stack[-1].shape[-1]]), 1)[:,0]
                logits = self.W_out(h_V_t) / temperature
                probs = jax.nn.softmax(logits-constant[None,:]*1e8+constant_bias[None,:]/temperature+bias_by_res_gathered/temperature, axis=-1)
                if pssm_bias_flag:
                    pssm_coef_gathered = jnp.take_along_axis(pssm_coef, t[:,None], 1)[:,0]
                    pssm_bias_gathered = jnp.take_along_axis(pssm_bias, t[:,None,None].tile([1,1,pssm_bias.shape[-1]]), 1)[:,0]
                    probs = (1-pssm_multi*pssm_coef_gathered[:,None])*probs + pssm_multi*pssm_coef_gathered[:,None]*pssm_bias_gathered
                if pssm_log_odds_flag:
                    pssm_log_odds_mask_gathered = jnp.take_along_axis(pssm_log_odds_mask, t[:,None, None].tile([1,1,pssm_log_odds_mask.shape[-1]]), 1)[:,0] #[B, 21]
                    probs_masked = probs*pssm_log_odds_mask_gathered
                    probs_masked += probs * 0.001
                    probs = probs_masked/jnp.sum(probs_masked, axis=-1, keepdims=True) #[B, 21]
                if omit_AA_mask_flag:
                    omit_AA_mask_gathered = jnp.take_along_axis(omit_AA_mask, t[:,None, None].tile([1,1,omit_AA_mask.shape[-1]]), 1)[:,0] #[B, 21]
                    probs_masked = probs*(1.0-omit_AA_mask_gathered)
                    probs = probs_masked/jnp.sum(probs_masked, axis=-1, keepdims=True) #[B, 21]
                used_key = jax.random.split(key, probs.shape[0])
                input = jnp.arange(probs.shape[1])[None].tile([probs.shape[0], 1])
                S_t = jax.vmap(lambda key, input, prob: jax.random.choice(key, input, p=prob),
                               in_axes=(0, 0, 0), out_axes=0)(used_key, input, probs)
                all_probs = scatter(all_probs, 1, t[:,None,None].tile([1,1,21]),
                                    jax.lax.convert_element_type((chain_mask_gathered[:,:,None,]*probs[:,None,:]),
                                                                 jnp.float32))
            S_true_gathered = jnp.take_along_axis(S_true, t[:,None], 1)
            S_t = jax.lax.convert_element_type((S_t*chain_mask_gathered+S_true_gathered*(1.0-chain_mask_gathered)),
                                               jnp.int64)
            temp1 = self.W_s(S_t)
            h_S = scatter(h_S, 1, t[:,None,None].tile([1,1,temp1.shape[-1]]), temp1)
            S = scatter(S, 1, t[:,None], S_t)
        output_dict = {"S": S, "probs": all_probs, "decoding_order": decoding_order}
        return output_dict


    def tied_sample(self, key, X, randn, S_true,
                    chain_mask, chain_encoding_all, residue_idx,
                    mask=None, temperature=1.0, omit_AAs_np=None,
                    bias_AAs_np=None, chain_M_pos=None, omit_AA_mask=None,
                    pssm_coef=None, pssm_bias=None, pssm_multi=None,
                    pssm_log_odds_flag=None, pssm_log_odds_mask=None,
                    pssm_bias_flag=None, tied_pos=None, tied_beta=None,
                    bias_by_res=None):

        # Prepare node and edge embeddings
        E, E_idx = self.features(X, mask, residue_idx, chain_encoding_all)
        h_V = jnp.zeros((E.shape[0], E.shape[1], E.shape[-1]))
        h_E = self.W_e(E)
        # Encoder is unmasked self-attention 
        mask_attend = gather_nodes(jnp.expand_dims(mask, -1),  E_idx).squeeze(-1)
        mask_attend = jnp.expand_dims(mask, -1) * mask_attend
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
        decoding_order = jnp.array(list(itertools.chain(*new_decoding_order)))[None,].tile([X.shape[0], 1])

        mask_size = E_idx.shape[1]
        permutation_matrix_reverse = jax.lax.convert_element_type(jax.nn.one_hot(decoding_order, num_classes=mask_size),
                                                                  jnp.float32)
        order_mask_backward = jnp.einsum('ij, biq, bjp->bqp',
                                         (1-jnp.triu(jnp.ones([mask_size,mask_size]))),
                                         permutation_matrix_reverse,
                                         permutation_matrix_reverse)
        mask_attend = jnp.expand_dims(jnp.take_along_axis(order_mask_backward, E_idx, 2), -1)
        mask_1D = mask.reshape([mask.shape[0], mask.shape[1], 1, 1])
        mask_bw = mask_1D * mask_attend
        mask_fw = mask_1D * (1. - mask_attend)

        N_batch, N_nodes = X.shape[0], X.shape[1]
        log_probs = jnp.zeros((N_batch, N_nodes, 21))
        all_probs = jnp.zeros((N_batch, N_nodes, 21), dtype=jnp.float32)
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
                input = jnp.arange(probs.shape[1])[None].tile([probs.shape[0], 1])
                S_t_repeat = jax.vmap(lambda key, input, prob: jax.random.choice(key, input, p=prob),
                                      in_axes=(0, 0, 0), out_axes=0)(used_key, input, probs)

                for t in t_list:
                    h_S = h_S.at[:,t,:].set(self.W_s(S_t_repeat))
                    S = S.at[:,t].set(S_t_repeat)
                    all_probs = all_probs.at[:,t,:].set(jax.lax.convert_element_type(probs, jnp.float32))
        output_dict = {"S": S, "probs": all_probs, "decoding_order": decoding_order}
        return output_dict
