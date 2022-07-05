from colabdesign.af.model import mk_af_model
from colabdesign.tr.model import mk_tr_model
class mk_af_tr_model:
  def __init__(self, protocol="fixbb", use_templates=False):
    assert protocol in ["fixbb","partial","hallucination","binder"]
    self.af = mk_afdesign_model(protocol=protocol, use_templates=use_templates)
    
    if protocol == "binder":
      def _prep_inputs(pdb_filename, chain, binder_len=50, binder_chain=None, **kwargs):
        self.af.prep_inputs(pdb_filename=pdb_filename, chain=chain,
                            binder_len=binder_len, binder_chain=binder_chain, **kwargs)
        if binder_chain is None:
          self.tr = mk_trdesign_model(protocol="hallucination")
          self.tr.prep_inputs(length=binder_len)
        else:
          self.tr = mk_trdesign_model(protocol="fixbb")
          self.tr.prep_inputs(pdb_filename=pdb_filename, chain=binder_chain)
    else:
      self.tr = mk_trdesign_model(protocol=protocol)

    if protocol == "fixbb":
      def _prep_inputs(pdb_filename, chain):
        self.af.prep_inputs(pdb_filename=pdb_filename, chain=chain)
        self.tr.prep_inputs(pdb_filename=pdb_filename, chain=chain)

    if protocol == "partial":
      def _prep_inputs(pdb_filename, chain, pos=None, length=None, fix_seq=False, use_sidechains=False):
        if use_sidechains: fix_seq = True
        flags = dict(pdb_filename=pdb_filename, chain=chain, pos=pos,
                     length=length, fix_seq=fix_seq)
        self.af.prep_inputs(**flags, use_sidechains=use_sidechains)
        self.tr.prep_inputs(**flags)    

    if protocol == "hallucintion":
      def _prep_inputs(length=None):
        self.af.prep_inputs(length=length)
        self.tr.prep_inputs(length=length)

    self.prep_inputs = _prep_inputs

  def set_opt(self,*args,**kwargs):
    self.af.set_opt(*args,**kwargs)
    self.tr.set_opt(*args,**kwargs)

  def predesign(self, pre_iter=100, logits=False, seed=None, verbose=10):
    self.tr.restart(seed=seed)
    self.tr.design(pre_iter, save_best=True, verbose=verbose)    
    seq = self.tr._best_aux["seq"]
    seq = seq["input"] if logits else seq["pseudo"]
    self.af.restart(seed=seed, seq=seq)