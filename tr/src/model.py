import jax.numpy as jnp
import jax
import numpy as np

def TrRosetta():  
  def pseudo_mrf(seq):
    '''single sequence'''
    L,A = seq.shape[0],21
    x_i = seq
    f_i = jnp.pad(seq,[[0,0],[0,1]])
    h_i = jnp.zeros((L,1))

    # 1D features
    feat_1D = jnp.concatenate([x_i,f_i,h_i],-1)
    feat_1D_tile_A = jnp.tile(feat_1D[:,None,:], [1,L,1])
    feat_1D_tile_B = jnp.tile(feat_1D[None,:,:], [L,1,1])

    # 2D features
    ic = 0.4 * jnp.eye(L*A)
    ic = jnp.reshape(ic,(L,A,L,A))
    ic = jnp.transpose(ic,(0,2,1,3))
    ic = jnp.reshape(ic,(L,L,A*A))
    i0 = jnp.zeros([L,L,1])
    feat_2D = jnp.concatenate([ic,i0],-1)

    feat = jnp.concatenate([feat_1D_tile_A, feat_1D_tile_B, feat_2D],-1)
    return feat[None]

  def instance_norm(x,p):
    mu = x.mean((1,2),keepdims=True)
    var = x.var((1,2),keepdims=True)
    inv = jax.lax.rsqrt(var + 1e-6) * p ["scale"]
    return x * inv + p["offset"] - mu * inv

  def conv_2D(x,p,d):
    c = dict(window_strides=(1,1), rhs_dilation=(d,d), padding="SAME")
    x = x.transpose([0,3,1,2])                # (batch, channels, row, col)
    f = p["filters"].transpose([3,2,0,1])     # (out_dims,in_dims,win,win)
    x = jax.lax.conv_general_dilated(x,f,**c)
    x = x.transpose([0,2,3,1])                # (batch, row, col, filters)
    return x + p["bias"]

  def dense(x,p):
    return x @ p["filters"] + p["bias"]
      
  def block(x,p,d):
    p = [jax.tree_map(lambda x:x[k],p) for k in [0,1]]
    y = conv_2D(x,p[0],d)
    y = instance_norm(y,p[0])
    y = jax.nn.elu(y)
    y = conv_2D(y,p[1],d)
    y = instance_norm(y,p[1])
    return jax.nn.elu(x + y)
  
  def resnet(x,p):
    for n,d in enumerate([1,2,4,8,16]):
      p_sub = jax.tree_map(lambda x:x[n],p)
      x = block(x,p_sub,d)
    return x, None

  def model(seq, params):
    x = pseudo_mrf(seq)
    x = dense(x,params["dense"])
    x = instance_norm(x,params["dense"])
    x = jax.nn.elu(x)
    x = jax.lax.scan(resnet,x,params["resnet"])[0]
    x = block(x,params["block"],1)
    
    o = {}
    o["theta"] = dense(x,params["theta"])
    o["phi"] = dense(x,params["phi"])
    x = (x + x.transpose(0,2,1,3))/2
    o["dist"] = dense(x,params["dist"])
    o["omega"] = dense(x,params["omega"])
    return jax.tree_map(lambda x:x[0], o)

  return model

def get_model_params(npy):
  '''parse TrRosetta params'''
  xaa = np.load(npy,allow_pickle=True)
  p,p_ = [],[[],[],[],[]]
  m = -1
  for n,w in enumerate(xaa):
    if w.shape == (3,3,64,64): m = 0
    if m > -1 and m < 4:
      p_[m].append(w)
      m += 1
    else:
      p.append(np.squeeze(w))
      m = -1
  p_ = [np.array(k) for k in p_]
  p0 = [k[:-2].reshape(12,5,2,*k.shape[1:]) for k in p_]
  p1 = [k[-2:].reshape(2,*k.shape[1:]) for k in p_]
  labels = ["filters","bias","offset","scale"]
  return {"dense": dict(zip(labels,p)),
          "resnet":dict(zip(labels,p0)),
          "block": dict(zip(labels,p1)),
          "theta": {"filters":p[4],"bias":p[5]},
          "phi":   {"filters":p[6],"bias":p[7]},
          "dist":  {"filters":p[8],"bias":p[9]},
          "omega": {"filters":p[12],"bias":p[13]}}
