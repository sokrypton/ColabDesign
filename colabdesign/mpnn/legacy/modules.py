import functools
import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np
import itertools
import joblib

from colabdesign.shared.prng import SafeKey
from .utils import gather_edges, gather_nodes, cat_neighbors_nodes, scatter, get_ar_mask
from .sample import mpnn_sample

Gelu = functools.partial(jax.nn.gelu, approximate=False)

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
    E = self.linear(d_onehot)
    return E

class RunModel:
  def __init__(self, config) -> None:
    self.config = config

    def _forward_score(inputs):
      model = ProteinMPNN(**self.config)
      return model(**inputs)
    self.score = jax.jit(hk.transform(_forward_score).apply)
    self.init_score = jax.jit(hk.transform(_forward_score).init)

    def _forward_sample(inputs):
      model = ProteinMPNN(**self.config)
      return model.sample(**inputs)
    self.sample = jax.jit(hk.transform(_forward_sample).apply)
    self.init_sample = jax.jit(hk.transform(_forward_sample).init)

    def _forward_tsample(inputs):
      model = ProteinMPNN(**self.config)
      return model.tied_sample(**inputs)
    self.tied_sample = jax.jit(hk.transform(_forward_tsample).apply)
    self.init_tsample = jax.jit(hk.transform(_forward_tsample).init)

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

  def _get_edge_idx(self, X, mask, eps=1E-6):
    ''' get edge index
    input: mask.shape = (...,L), X.shape = (...,L,3)
    return: (...,L,k)
    '''
    mask_2D = mask[...,None,:] * mask[...,:,None] 
    dX = X[...,None,:,:] - X[...,:,None,:]
    D = jnp.sqrt(jnp.square(dX).sum(-1) + eps)
    D_masked = jnp.where(mask_2D,D,D.max(-1,keepdims=True))
    k = min(self.top_k, X.shape[-2])
    _, E_idx = jax.lax.approx_min_k(D_masked, k, reduction_dimension=-1)
    return E_idx

  def _rbf(self, D):
    ''' radial basis function (RBF)
    input: (...,L,k)
    output: (...,L,k,?)
    '''
    D_min, D_max, D_count = 2., 22., self.num_rbf
    D_mu = jnp.linspace(D_min, D_max, D_count)
    D_sigma = (D_max - D_min) / D_count    
    return jnp.exp(-((D[...,None] - D_mu) / D_sigma)**2)

  def _get_rbf(self, A, B, E_idx):
    D = jnp.sqrt(jnp.square(A[...,:,None,:] - B[...,None,:,:]).sum(-1) + 1e-6)
    D_neighbors = gather_edges(D[...,None], E_idx)[...,0] #[...,L,K]
    return self._rbf(D_neighbors)

  def __call__(self, X, mask, residue_idx, chain_idx, offset=None):
    if self.augment_eps > 0:
      self.safe_key, use_key = self.safe_key.split()
      X = X + self.augment_eps * jax.random.normal(use_key, X.shape)
    
    ##########################
    # get atoms
    ##########################
    # N,Ca,C,O,Cb
    Y = X.transpose((2,0,1,3))
    if Y.shape[0] == 4:
      # add Cb
      b,c = (Y[1]-Y[0]),(Y[2]-Y[1])
      Cb = -0.58273431*jnp.cross(b,c) + 0.56802827*b - 0.54067466*c + Y[1]
      Y = jnp.concatenate([Y,Cb[None]],0)

    ##########################
    # gather edge features
    ##########################
    # get edge indices (based on ca-ca distances)
    E_idx = self._get_edge_idx(Y[1], mask)

    # rbf encode distances between atoms
    edges = jnp.array([[1,1],[0,0],[2,2],[3,3],[4,4],
                       [1,0],[1,2],[1,3],[1,4],[0,2],
                       [0,3],[0,4],[4,2],[4,3],[3,2],
                       [0,1],[2,1],[3,1],[4,1],[2,0],
                       [3,0],[4,0],[2,4],[3,4],[2,3]])
    RBF_all = jax.vmap(lambda x:self._get_rbf(Y[x[0]],Y[x[1]],E_idx))(edges)
    RBF_all = RBF_all.transpose((1,2,3,0,4))
    RBF_all = RBF_all.reshape(RBF_all.shape[:-2]+(-1,))

    ##########################
    # position embedding
    ##########################
    # residue index offset
    if offset is None:
      offset = (residue_idx[...,:,None] - residue_idx[...,None,:])
    offset = gather_edges(offset[...,None], E_idx)[...,0] #[B, L, K]

    # chain index offset
    d_chains = (chain_idx[...,:,None] == chain_idx[...,None,:]).astype(int)
    E_chains = gather_edges(d_chains[...,None], E_idx)[...,0]
    E_positional = self.embeddings(offset, E_chains)

    ##########################
    # define edges
    ##########################
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
    if jnp.issubdtype(arr.dtype, jnp.integer):
      one_hot = jax.nn.one_hot(arr, self.vocab_size)
    else:
      one_hot = arr
    return jnp.tensordot(one_hot, self.embeddings, 1)

class ProteinMPNN(hk.Module, mpnn_sample):
  def __init__(self, num_letters,
         node_features, edge_features, hidden_dim,
         num_encoder_layers=3, num_decoder_layers=3,
         vocab=21, k_neighbors=64,
         augment_eps=0.05, dropout=0.1):
    super(ProteinMPNN, self).__init__()

    # Hyperparameters
    self.node_features = node_features
    self.edge_features = edge_features
    self.hidden_dim = hidden_dim

    # Featurization layers
    self.features = ProteinFeatures(edge_features,
                    node_features,
                    top_k=k_neighbors,
                    augment_eps=augment_eps)

    self.W_e = hk.Linear(hidden_dim, with_bias=True, name='W_e')
    self.W_s = EmbedToken(vocab_size=vocab, embed_dim=hidden_dim)

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

  def __call__(self, X, mask, residue_idx, chain_idx,
               S=None, chain_M=None, randn=None,
               ar_mask=None, decoding_order=None, offset=None):
    """ Graph-conditioned sequence model """
    # Prepare node and edge embeddings
    E, E_idx = self.features(X, mask, residue_idx, chain_idx, offset=offset)
    h_V = jnp.zeros((E.shape[0], E.shape[1], E.shape[-1]))
    h_E = self.W_e(E)

    # Encoder is unmasked self-attention
    mask_attend = gather_nodes(mask[...,None],E_idx)[...,0]
    mask_attend = mask[...,None] * mask_attend
    for layer in self.encoder_layers:
      h_V, h_E = layer(h_V, h_E, E_idx, mask, mask_attend)

    # Build encoder embeddings
    h_EX_encoder = cat_neighbors_nodes(jnp.zeros_like(h_V), h_E, E_idx)
    h_EXV_encoder = cat_neighbors_nodes(h_V, h_EX_encoder, E_idx)

    if S is None:  
      ##########################################
      # unconditional_probs
      ##########################################

      # make an autogressive mask
      ar_mask = jnp.zeros([X.shape[0], X.shape[1], X.shape[1]])

      mask_attend = jnp.take_along_axis(ar_mask, E_idx, 2)[...,None]
      mask_1D = mask.reshape([mask.shape[0], mask.shape[1], 1, 1])
      mask_bw = mask_1D * mask_attend
      mask_fw = mask_1D * (1. - mask_attend)

      h_EXV_encoder_fw = mask_fw * h_EXV_encoder
      for layer in self.decoder_layers:
        h_V = layer(h_V, h_EXV_encoder_fw, mask)

    else:
      ##########################################
      # conditional_probs
      ##########################################

      # Concatenate sequence embeddings for autoregressive decoder
      h_S = self.W_s(S)
      h_ES = cat_neighbors_nodes(h_S, h_E, E_idx)

      if ar_mask is None:
        if decoding_order is None:
          # update chain_M to include missing regions
          chain_M = chain_M * mask
          #[numbers will be smaller for places where chain_M = 0.0 and higher for places where chain_M = 1.0]
          decoding_order = jnp.argsort((chain_M+0.0001)*(jnp.abs(randn)))
        
        # make an autogressive mask
        ar_mask = get_ar_mask(decoding_order)
              
      mask_attend = jnp.take_along_axis(ar_mask, E_idx, 2)[...,None]
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