from colabdesign.shared.utils import update_dict
class _model:
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