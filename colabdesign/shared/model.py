from colabdesign.shared.utils import update_dict
from colabdesign.shared.prep import rewire

class design_model:
  def set_opt(self, *args, **kwargs):
    '''set [opt]ions'''
    if kwargs.pop("set_defaults", False):
      update_dict(self._opt, *args, **kwargs)
    update_dict(self.opt, *args, **kwargs)

  def set_weights(self, *args, **kwargs):
    '''set weights'''
    if kwargs.pop("set_defaults", False):
      update_dict(self._opt["weights"], *args, **kwargs)
    update_dict(self.opt["weights"], *args, **kwargs)

  def rewire(self, order=None, offset=0, loops=0, set_defaults=True):
    if "pos" in self.opt:
      self.opt["pos"] = rewire(length=self._pos_info["length"], order=order,
                               offset=offset, loops=loops)
      if set_defaults: self._opt["pos"] = self.opt["pos"]