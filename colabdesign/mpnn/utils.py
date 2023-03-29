import jax.numpy as jnp
import jax

def gather_nodes(nodes, neighbor_idx):
  # Features [B,N,C] at Neighbor indices [B,N,K] => [B,N,K,C]
  # Flatten and expand indices per batch [B,N,K] => [B,NK] => [B,NK,C]
  neighbors_flat = neighbor_idx.reshape([neighbor_idx[None].shape[0], -1])
  neighbors_flat = jnp.tile(jnp.expand_dims(neighbors_flat, -1),[1, 1, nodes[None].shape[2]])
  # Gather and re-pack
  neighbor_features = jnp.take_along_axis(nodes[None], neighbors_flat, 1)
  neighbor_features = neighbor_features.reshape(list(neighbor_idx[None].shape[:3]) + [-1])
  return neighbor_features[0]

def cat_neighbors_nodes(h_nodes, h_neighbors, E_idx):
  h_nodes = gather_nodes(h_nodes, E_idx)[None]
  h_nn = jnp.concatenate([h_neighbors[None], h_nodes], -1)
  return h_nn[0]

def get_ar_mask(order):
  '''compute autoregressive mask, given order of positions'''
  order = order.flatten()
  L = order.shape[-1]
  tri = jnp.tri(L, k=-1)
  idx = order.argsort()
  ar_mask = tri[idx,:][:,idx]
  return ar_mask