import numpy as np
import itertools
import jax.numpy as jnp
import jax
import time

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

def get_ar_mask(order):
  '''compute autoregressive mask, given order of positions'''
  L = order.shape[-1]
  oh_order = jax.nn.one_hot(order, L)
  tri = jnp.tri(L, k=-1)
  return jnp.einsum('ij,...iq,...jp->...qp', tri, oh_order, oh_order)