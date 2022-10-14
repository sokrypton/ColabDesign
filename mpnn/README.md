# ProteinMPNN in jax!
**WARNING** This code is work-in-progress!

<a href="https://colab.research.google.com/github/sokrypton/ColabDesign/blob/v1.1.0/mpnn/examples/proteinmpnn_in_jax.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

#### install
```bash
pip -q install git+https://github.com/sokrypton/ColabDesign.git@v1.1.0
```
#### run
```python
from colabdesign.mpnn import mk_mpnn_model
mpnn_model = mk_mpnn_model()
mpnn_model.prep_inputs(pdb_filename=???, chain=???)
mpnn_model.sample_parallel(batch=128)
```
### Contributors:
- Sergey Ovchinnikov [@sokrypton](https://github.com/sokrypton)
- Shihao Feng [@JeffSHF](https://github.com/JeffSHF)
- Simon Kozlov [@sim0nsays](https://github.com/sim0nsays)
- Justas Dauparas [@dauparas](https://github.com/dauparas) - [original pytorch code](https://github.com/dauparas/ProteinMPNN)
