from colabdesign.shared.utils import update_dict
class _model:
  def set_opt(self, *args, **kwargs):
    '''set [opt]ions'''
    update_dict(self.opt, *args, **kwargs)
    if kwargs.pop("set_default", False):
      update_dict(self._opt, *args, **kwargs)

  def set_weights(self, *args, **kwargs):
    '''set weights'''
    update_dict(self.opt["weights"], *args, **kwargs)
    if kwargs.pop("set_default", False):
      update_dict(self._opt["weights"], *args, **kwargs)