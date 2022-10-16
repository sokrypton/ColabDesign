from colabdesign.af.model import mk_af_model
from colabdesign.tr.model import mk_tr_model

class mk_af_tr_model:
  def __init__(self, protocol="fixbb", use_templates=False,
               recycle_mode="last", num_recycles=0):
    assert protocol in ["fixbb","partial","hallucination","binder"]
    self.af = mk_af_model(protocol=protocol, use_templates=use_templates,
                          recycle_mode=recycle_mode, num_recycles=num_recycles)
    
    if protocol == "binder":
      def _prep_inputs(pdb_filename, chain, binder_len=50, binder_chain=None,
                       ignore_missing=True, **kwargs):
        self.af.prep_inputs(pdb_filename=pdb_filename, chain=chain,
                            binder_len=binder_len, binder_chain=binder_chain,
                            ignore_missing=ignore_missing, **kwargs)
        flags = dict(ignore_missing=ignore_missing)
        if binder_chain is None:
          self.tr = mk_tr_model(protocol="hallucination")
          self.tr.prep_inputs(length=binder_len, **flags)
        else:
          self.tr = mk_tr_model(protocol="fixbb")
          self.tr.prep_inputs(pdb_filename=pdb_filename, chain=binder_chain, **flags)
    else:
      self.tr = mk_tr_model(protocol=protocol)

    if protocol == "fixbb":
      def _prep_inputs(pdb_filename, chain, fix_pos=None, 
                       ignore_missing=True, **kwargs):
        flags = dict(pdb_filename=pdb_filename, chain=chain,
                     fix_pos=fix_pos, ignore_missing=ignore_missing)
        self.af.prep_inputs(**flags, **kwargs)
        self.tr.prep_inputs(**flags, chain=chain)

    if protocol == "partial":
      def _prep_inputs(pdb_filename, chain, pos=None, length=None,
                       fix_pos=None, use_sidechains=False, atoms_to_exclude=None, 
                       ignore_missing=True, **kwargs):
        if use_sidechains: fix_seq = True
        flags = dict(pdb_filename=pdb_filename, chain=chain, 
                     length=length, pos=pos, fix_pos=fix_pos,
                     ignore_missing=ignore_missing)
        af_a2e = kwargs.pop("af_atoms_to_exclude",atoms_to_exclude)
        tr_a2e = kwargs.pop("tr_atoms_to_exclude",atoms_to_exclude)
        self.af.prep_inputs(**flags, use_sidechains=use_sidechains, atoms_to_exclude=af_a2e, **kwargs)
        self.tr.prep_inputs(**flags, atoms_to_exclude=tr_a2e)
      
      def _rewire(order=None, offset=0, loops=0):
        self.af.rewire(order=order, offset=offset, loops=loops)
        self.tr.rewire(order=order, offset=offset, loops=loops)
      
      self.rewire = _rewire

    if protocol == "hallucintion":
      def _prep_inputs(length=None, **kwargs):
        self.af.prep_inputs(length=length, **kwargs)
        self.tr.prep_inputs(length=length)

    self.prep_inputs = _prep_inputs

  def set_opt(self,*args,**kwargs):
    self.af.set_opt(*args,**kwargs)
    self.tr.set_opt(*args,**kwargs)

  def joint_design(self, iters=100, tr_weight=1.0, tr_seed=None, **kwargs):
    self.af.design(iters, callback=self.tr.af_callback(weight=tr_weight, seed=tr_seed), **kwargs)
