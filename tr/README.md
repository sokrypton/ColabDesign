# TrDesign in JAX!
More info to come! Currently only `fixbb` and `hallucination` protocols are supported.

Download weights
```bash
if [ ! -d models ]; then
  wget -qnc https://files.ipd.uw.edu/krypton/TrRosetta/models.zip
  wget -qnc https://files.ipd.uw.edu/krypton/TrRosetta/bkgr_models.zip
  unzip -qqo models.zip
  unzip -qqo bkgr_models.zip
fi
```

Example
```python
from af import mk_afdesign_model, clear_mem
from tr import mk_trdesign_model

clear_mem()
af_model = mk_afdesign_model(protocol="fixbb")
af_model.prep_inputs(get_pdb("1TEN"))

tr_model = mk_trdesign_model(protocol="fixbb")
tr_model.prep_inputs(get_pdb("1TEN"))

af_model.restart()
af_model.design_3stage(callback=tr_model.af_callback(weight=1.0))
```

