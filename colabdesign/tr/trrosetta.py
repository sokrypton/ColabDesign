import jax.numpy as jnp
import jax
import numpy as np

def TrRosetta(bkg_model=False):  
  
  def pseudo_mrf(inputs, prf=None):
    '''single sequence'''
    seq,prf = inputs["seq"],inputs["prf"]
    L,A = seq.shape[0],21
    if prf.shape[1] == 20:
      prf = jnp.pad(prf,[[0,0],[0,1]])
    
    # 1D features
    x_1D = jnp.concatenate([seq, prf],-1)
    x_1D = jnp.pad(x_1D,[[0,0],[0,1]])
    x_1D = jnp.repeat(x_1D[None],L,0)

    # 2D features
    x_2D = jnp.diag(jnp.full(L*A,0.4))
    x_2D = x_2D.reshape(L,A,L,A).swapaxes(1,2).reshape(L,L,-1)
    x_2D = jnp.pad(x_2D,[[0,0],[0,0],[0,1]])
    return jnp.concatenate([x_1D.swapaxes(0,1), x_1D, x_2D],-1)
  
  # layers
  def instance_norm(x, params):
    mu = x.mean((0,1),keepdims=True)
    var = x.var((0,1),keepdims=True)
    inv = jax.lax.rsqrt(var + 1e-6) * params["scale"]
    return x * inv + params["offset"] - mu * inv

  def conv_2D(x, params, dilation=1, stride=1, padding="SAME"):
    flags = dict(window_strides=(stride,stride),
                 rhs_dilation=(dilation,dilation),
                 padding=padding)
    x = x.transpose([2,0,1])
    f = params["filters"].transpose([3,2,0,1])
    x = jax.lax.conv_general_dilated(x[None], f, **flags)[0]
    x = x.transpose([1,2,0])
    return x + params["bias"]

  def dense(x, params):
    return x @ params["filters"] + params["bias"]

  def dropout(x, key, rate):
    keep_rate = 1.0 - rate
    keep = jax.random.bernoulli(key, keep_rate, shape=x.shape)
    return keep * x / keep_rate
  
  # meta layers
  def encoder(x, params):
    x = dense(x, params)
    x = instance_norm(x, params)
    return jax.nn.elu(x)
      
  def block(x, params, dilation, key, rate=0.15):
    y = x
    for n in [0,1]:
      if n == 1: y = dropout(y, key, rate)
      p = jax.tree_map(lambda x:x[n], params)
      y = conv_2D(y, p, dilation)
      y = instance_norm(y, p)
      y = jax.nn.elu(y if n == 0 else (x+y))
    return y  
  
  def resnet(x, params, key, rate=0.15):
    def body(prev, sub_params):
      (x,key) = prev
      for n, dilation in enumerate([1,2,4,8,16]):
        key, sub_key = jax.random.split(key)
        p = jax.tree_map(lambda x:x[n], sub_params)
        x = block(x, p, dilation, sub_key, rate)
      return (x,key), None
    return jax.lax.scan(body,(x,key),params)[0][0]
  
  def heads(x, params):
    o = {k:dense(x,params[k]) for k in ["theta","phi"]}
    x = (x + x.swapaxes(0,1)) / 2  
    o.update({k:dense(x,params[k]) for k in ["dist","bb","omega"]})
    return o    
  
  def trunk(x, params, key, rate=0.15):
    key, sub_key = jax.random.split(key)
    x = encoder(x, params["encoder"])
    x = resnet(x, params["resnet"], sub_key, rate)
    x = block(x, params["block"], 1, key, rate)  
    return heads(x, params)
  
  # decide which model to use
  if bkg_model:
    def model(params, key, length=100):
      key, sub_key = jax.random.split(key)
      x = jax.random.normal(sub_key, (length, length, 64))
      return trunk(x, params, key, 0.0)
    return jax.jit(model, static_argnums=2)
  else:
    def model(inputs, params, key, rate=0.15):
      x = pseudo_mrf(inputs)
      return trunk(x, params, key, rate)
    return jax.jit(model)

def get_model_params(npy):
  '''parse TrRosetta params into dictionary'''
  xaa = np.load(npy,allow_pickle=True).tolist()
  layers = ["encoder","resnet","block","theta","phi","dist","bb","omega"]
  num = np.array([4,0,8,2,2,2,2,2])
  num[1] = len(xaa) - num.sum()
  idx = np.cumsum(num) - num
  def split(params):
    labels = ["filters","bias","offset","scale"]
    steps = min(len(params),len(labels))
    return {labels[n]:np.squeeze(params[n::steps]) for n in range(steps)}
  params = {k:split(xaa[i:i+n]) for k,i,n in zip(layers,idx,num)}
  params["resnet"] = jax.tree_map(lambda x:x.reshape(-1,5,2,*x.shape[1:]), params["resnet"])
  return params 