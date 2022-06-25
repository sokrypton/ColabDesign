def update_dict(d, u=None):
  if u is not None:
    for k, v in u.items():
      if v is not None:
        if isinstance(v, dict):
          update_dict(d.get(k,{}), v)
        else:
          d.update({k:v})
