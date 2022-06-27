import jax
def update_dict(D, *args, **kwargs):
  def set_dict(d, x):
    for k,v in x.items():
      if v is not None:
        if k in d:
          if isinstance(v,dict):
            set_dict(d[k], x[k])
          elif isinstance(d[k],dict):
            print(f"ERROR: '{k}' is a dictionary")
          elif isinstance(d[k],(int,float,bool,str)):
            d[k] = type(d[k])(v)
          else:
            d[k] = v
        else:
          print(f"ERROR: '{k}' not found in {list(d.keys())}")  
  for a in args:
    if isinstance(a, dict): set_dict(D, a)    
  set_dict(D, kwargs)
