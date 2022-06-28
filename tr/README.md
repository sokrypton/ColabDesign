# TrDesign in JAX!
Work in Progress...
For original version implemented in keras see: https://github.com/gjoni/trDesign/tree/master/02-GD

### Google Colab
<a href="https://colab.research.google.com/github/sokrypton/ColabDesign/blob/main/af/design.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

### install
```bash
pip install git+https://github.com/sokrypton/ColabDesign.git

# download weights
if [ ! -d params/tr ]; then
  mkdir -p params/tr
  wget https://files.ipd.uw.edu/krypton/TrRosetta/models.zip
  wget https://files.ipd.uw.edu/krypton/TrRosetta/bkgr_models.zip
  unzip models.zip -d params/tr/
  unzip bkgr_models.zip -d params/tr/
fi
```

### example
```python
from colabdesign import clear_mem, mk_trdesign_model

clear_mem()
tr_model = mk_trdesign_model(protocol="fixbb")
tr_model.prep_inputs(get_pdb("6MRR"), chain="A")
tr_model.design(100, verbose=10)
tr_model.plot()
print(tr_model.get_loss())
print(tr_model.get_seq())
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

