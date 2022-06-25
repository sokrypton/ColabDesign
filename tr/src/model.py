import jax.numpy as jnp
import jax
import numpy as np

def TrRosetta(mode="Tr"):  
  
  assert mode in ["Tr","Tr_bkg"]
  
  def pseudo_mrf(seq):
    '''single sequence'''    
    L,A = seq.shape[0],21
    prf = jnp.pad(seq,[[0,0],[0,1]])
    
    # 1D features
    x_1D = jnp.concatenate([seq, prf],-1)
    x_1D = jnp.pad(x_1D,[[0,0],[0,1]])
    x_1D = jnp.repeat(x_1D[None],L,0)

    # 2D features
    x_2D = jnp.diag(jnp.full(L*A,0.4))
    x_2D = x_2D.reshape(L,A,L,A).swapaxes(1,2).reshape(L,L,-1)
    x_2D = jnp.pad(x_2D,[[0,0],[0,0],[0,1]])
    return jnp.concatenate([x_1D.swapaxes(0,1), x_1D, x_2D],-1)
  
  def instance_norm(x, p):
    mu = x.mean((0,1),keepdims=True)
    var = x.var((0,1),keepdims=True)
    inv = jax.lax.rsqrt(var + 1e-6) * p ["scale"]
    return x * inv + p["offset"] - mu * inv

  def conv_2D(x, p, dilation=1, stride=1, padding="SAME"):
    flags = dict(window_strides=(stride,stride),
                 rhs_dilation=(dilation,dilation),
                 padding=padding)
    x = x.transpose([2,0,1])
    f = p["filters"].transpose([3,2,0,1])
    x = jax.lax.conv_general_dilated(x[None], f, **flags)[0]
    x = x.transpose([1,2,0])
    return x + p["bias"]

  def dense(x, p):
    return x @ p["filters"] + p["bias"]
  
  def encoder(x, p):
    x = dense(x,p)
    x = instance_norm(x,p)
    return jax.nn.elu(x)
      
  def block(x, p, dilation=1):
    y = x
    for n in [0,1]:
      q = jax.tree_map(lambda x:x[n],p)
      y = conv_2D(y,q,dilation)
      y = instance_norm(y,q)
      y = jax.nn.elu(y if n == 0 else (x+y))
    return y  
  
  def resnet(x, p):
    def body(y, p_sub):
      for n,dilation in enumerate([1,2,4,8,16]):
        q = jax.tree_map(lambda x:x[n],p_sub)
        y = block(y,q,dilation)
      return y, None
    return jax.lax.scan(body,x,p)[0]
  
  def heads(x, p):
    o = {k:dense(x,p[k]) for k in ["theta","phi"]}
    x = (x + x.swapaxes(0,1)) / 2  
    o.update({k:dense(x,p[k]) for k in ["dist","bb","omega"]})
    return o    
  
  if mode == "Tr":
    def model(seq, model_params):
      x = pseudo_mrf(seq)
      x = encoder(x, model_params["encoder"])
      x = resnet(x, model_params["resnet"])
      x = block(x, model_params["block"])  
      x = heads(x, model_params)
    return model
  
  if mode == "Tr_bkg":
    def model(length, model_params, seed):
      x = jax.random.normal(seed, (length, length, 64))
      x = encoder(x, model_params["encoder"])
      x = resnet(x, model_params["resnet"])
      x = block(x, model_params["block"])  
      x = heads(x, model_params)
    return model

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
