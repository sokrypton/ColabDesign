# ProteinMPNN in jax!
**WARNING** This code is work-in-progress!

<a href="https://colab.research.google.com/github/sokrypton/ColabDesign/blob/v1.1.1/mpnn/examples/proteinmpnn_in_jax.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

#### install
```bash
pip -q install git+https://github.com/sokrypton/ColabDesign.git@v1.1.1
```
#### run
```python
from colabdesign.mpnn import mk_mpnn_model
mpnn_model = mk_mpnn_model()
mpnn_model.prep_inputs(pdb_filename="tmp.pdb")
samples = mpnn_model.sample_parallel()
```
# FAQ
#### What are all the available functions?
- `mpnn_model.sample()` - sample one sequence
- `mpnn_model.sample(temperature=0.1)` - control sampling temperature
- `mpnn_model.sample(decoding_order=np.array([0,1,2,3,4,5]))` specify the order of autoregressive sampling
- `mpnn_model.sample(decoding_order=np.array([[0,3],[1,4],[2,5]]))` - specify order of "tied" autoregressive sampling
- `mpnn_model.sample_parallel(batch=128)` - sample 128 sequences in parallel (all options above apply)
- `mpnn_model.score(seq="QWERTY")` - score one sequence
- `mpnn_model.get_unconditional_logits()` - get P(sequence | structure)
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
You can modify the bias matrix directly! The bias matrix is a (length, 21) matrix. Using large negative/positive values in the bias matrix is how we prevent certain amino acids from being sampled (rm_aa) and fix certain positions (fix_pos). For reference, the alphabet used: `ARNDCQEGHILKMFPSTWYV`.

For example, to add alanine bias to the first position, do:
```python
from colabdesign.mpnn.model import aa_order
mpnn_model.prep_inputs(pdb_filename="tmp.pdb")
mpnn_model._inputs["bias"][0,aa_order["A"]] = 1.0
```
For example, if you want to add a hydrophilic bias to all positions, you can do:
```python
for k in "DEHKNQRSTWY":
  mpnn_model._inputs["bias"][:,aa_order[k]] += 1.39
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
- Shihao Feng [@JeffSHF](https://github.com/JeffSHF)
- Sergey Ovchinnikov [@sokrypton](https://github.com/sokrypton)
- Simon Kozlov [@sim0nsays](https://github.com/sim0nsays)
- Justas Dauparas [@dauparas](https://github.com/dauparas) - [original pytorch code](https://github.com/dauparas/ProteinMPNN)
