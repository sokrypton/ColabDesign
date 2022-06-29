import jax
import numpy as jnp
import jax.numpy as jnp

def clear_mem():
  backend = jax.lib.xla_bridge.get_backend()
  for buf in backend.live_buffers(): buf.delete()

def update_dict(D, *args, **kwargs):
  '''robust function for updating dictionary'''
  def set_dict(d, x):
    for k,v in x.items():
      if v is not None:
        if k in d:
          if isinstance(v, dict):
            set_dict(d[k], x[k])
          elif isinstance(d[k],(np.ndarray,jnp.ndarray)):
            d[k] = np.asarray(v)
          else:
            d[k] = type(d[k])(v)
        else:
          print(f"ERROR: '{k}' not found in {list(d.keys())}")  
  while len(args) > 0 and isinstance(args[0],str):
    D,args = D[args[0]],args[1:]
  for a in args:
    if isinstance(a, dict): set_dict(D, a)
  set_dict(D, kwargs)