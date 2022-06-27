# TrDesign in JAX!
Work in Progress...
Currently only `fixbb` and `hallucination` protocols are supported as callbacks to AfDesign.
Experimental support for `partial` added.

For original version implemented in keras see: https://github.com/gjoni/trDesign/tree/master/02-GD

### download weights
```bash
%%bash
if [ ! -d models ]; then
  wget -qnc https://files.ipd.uw.edu/krypton/TrRosetta/models.zip
  wget -qnc https://files.ipd.uw.edu/krypton/TrRosetta/bkgr_models.zip
  unzip -qqo models.zip
  unzip -qqo bkgr_models.zip
fi
```

### example
combine AfDesign and TrDesign for fixed backbone design 
```python
from colabdesign import clear_mem, mk_afdesign_model, mk_trdesign_model

clear_mem()
af_model = mk_afdesign_model(protocol="fixbb")
af_model.prep_inputs(get_pdb("1TEN"))

tr_model = mk_trdesign_model(protocol="fixbb")
tr_model.prep_inputs(get_pdb("1TEN"))

af_model.restart()
af_model.design_3stage(callback=tr_model.af_callback())
```

