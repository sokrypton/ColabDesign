# TrDesign in JAX!
More info to come!
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
