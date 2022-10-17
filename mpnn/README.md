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
mpnn_model.prep_inputs(pdb_filename="tmp.pdb")
samples = mpnn_model.sample_parallel()
```
# FAQ
#### How do I specify which positions to fix, while leaving the rest to redesign?
```python
mpnn_model.prep_inputs(pdb_filename="tmp.pdb", fix_pos="1-10")
```
#### Can I invert the selection? So I can specify which positions to redesign?
```python
mpnn_model.prep_inputs(pdb_filename="tmp.pdb", fix_pos="1-10", inverse=True)
```
#### How about multichain inputs?
```python
mpnn_model.prep_inputs(pdb_filename="tmp.pdb", chain="A,B", fix_pos="A1-10,B5-20")
```
#### Can I fix an entire chain, for binder redesign?
```python
mpnn_model.prep_inputs(pdb_filename="tmp.pdb", chain="A,B", fix_pos="A")
```
#### Can I avoid certain amino acids?
```python
mpnn_model.prep_inputs(pdb_filename="tmp.pdb", rm_aa="C")
```
#### I want more control!
You can modify the bias matrix directly! The bias matrix is a (length, 20) matrix. Setting values to large negative will avoid these from being sampled (this is what rm_aa does) and setting values to large positions numbers will fix those positions (this is what fix_pos does). The alphabet for the 20 values is `ARNDCQEGHILKMFPSTWYV`.
```python
mpnn_model._inputs["bias"][:,0] = 1e8
```
#### How about tied sampling for homo-oligomeric complexes?
```python
mpnn_model.prep_inputs(pdb_filename="tmp.pdb", chain="A,B,C", homooligomeric=True)
```
# Advanced FAQ
#### How do I evaluate the sequences with AlphaFold?
```bash
mkdir params
curl -fsSL https://storage.googleapis.com/alphafold/alphafold_params_2022-03-02.tar | tar x -C params
```
```python
from colabdesign.af import mk_af_model
af_model = mk_af_model()
af_model.prep_inputs(pdb_filename="tmp.pdb")
for n,S in enumerate(samples["S"]):
  af_model.predict(seq=S.argmax(-1))
  af_model.save_current_pdb(f"{n}.pdb")
```

### Contributors:
- Sergey Ovchinnikov [@sokrypton](https://github.com/sokrypton)
- Shihao Feng [@JeffSHF](https://github.com/JeffSHF)
- Simon Kozlov [@sim0nsays](https://github.com/sim0nsays)
- Justas Dauparas [@dauparas](https://github.com/dauparas) - [original pytorch code](https://github.com/dauparas/ProteinMPNN)
