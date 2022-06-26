import jax

def update_dict(D, *args, **kwargs):
  def set_zero(x):
    for k,v in x.items():
      if isinstance(v,dict): set_zero(v)
      else: x[k] = type(x[k])(False)
  
  def set_dict(d, x):
    for k,v in x.items():
      if k in d:
        if isinstance(v,dict):
          set_dict(d[k], x[k])
        elif isinstance(d[k],dict):
          print(f"ERROR: '{k}' is a dictionary")
        elif type(d[k]) in [int,float,bool,str]:
          d[k] = type(d[k])(v)
        else:
          d[k] = v
      else:
        print(f"ERROR: '{k}' not found in {list(d.keys())}")
  
  if "zero_all" in kwargs and kwargs["zero_all"]:
    set_zero(D)
    kwargs.pop("zero_all")
  
  for a in args:
    if isinstance(a, dict): set_dict(D, a)    
  
  set_dict(D, kwargs)
