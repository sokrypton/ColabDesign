# AfDesign (v1.0.4)
### Google Colab
<a href="https://colab.research.google.com/github/sokrypton/ColabDesign/blob/main/af/design.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

# Updates
- **24Feb2022** - "Beta" branch started. Refactoring code to allow homooligomeric hallucination/design and averaging gradients across recycles (which is now the default).
Minor changes changes include renaming intra_pae/inter_con to pae/con and inter_pae/inter_con to i_pae/i_con for clarity.
- **28Feb2022** - We find backprop through structure module to be unstable, all functions have been updated to only use distogram by default. The definition of contact has changed to minimize entropy within distance cutoff.
- **02May2022** - The `design.py` code has been split up into multiple python files under `src/`
- **14May2022** - Adding support for partial hallucination (if you want to constrain one part and generate structure/sequence for rest).
- **19June2022** - "Beta" branch is now the "Main" branch. WARNING: Lots of default settings and weights were changed. [Click here](#i-was-getting-better-results-before-the-major-update-19june2022-how-do-i-revert-back-to-the-old-settings) for info on how to revert back to old settings. 
- **28June2022** - v1.0.1 - Major code reorganization/refactoring to add support for callbacks (to allow integration w/ other tools during design) and to avoid clashes with existing trrosetta/alphafold installations. (eg. `af → colabdesign`, `af.src → colabdesign.af` and `alphafold → colabdesign.af.alphafold`).
- **05July2022** - v1.0.2 - Major code cleanup, removing duplicate code. Adding support for custom loss functions.
- **11July2022** - v1.0.3 - Improved homo-oligomeric support. RMSD and dgram losses have been refactored to automatically save aligned coordinates. Multimeric coordinates now saved with chain identifiers.
- **23July2022** - v1.0.4 - Refactoring... Adding support for openfold weights. To enable set `mk_afdesign_model(...,use_openfold=True)`.

### setup
```bash
pip install git+https://github.com/sokrypton/ColabDesign.git

# download alphafold weights
mkdir params
curl -fsSL https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar | tar x -C params

# download openfold weights (optional)
for W in openfold_model_ptm_1 openfold_model_ptm_2 openfold_model_no_templ_ptm_1
do wget -qnc https://files.ipd.uw.edu/krypton/openfold/${W}.npz -P params; done
```
### import
```python
import numpy as np
from IPython.display import HTML
from colabdesign import mk_afdesign_model, clear_mem
```
### fixed backbone design
For a given protein backbone, generate/design a new sequence that AlphaFold thinks folds into that conformation
```python
model = mk_afdesign_model(protocol="fixbb")
model.prep_inputs(pdb_filename="1TEN.pdb", chain="A")
model.design_3stage()
```
### hallucination
For a given length, generate/hallucinate a protein sequence that AlphaFold thinks folds into a well structured 
protein (high plddt, low pae, many contacts).
```python
model = mk_afdesign_model(protocol="hallucination")
model.prep_inputs(length=100)
model.restart(seq_init="gumbel")
model.design(50, soft=True)
model.restart(seq_init=model.aux["seq"]["pseudo"], keep_history=True)
model.design_3stage(50,50,10)
```
### binder hallucination
For a given protein target and protein binder length, generate/hallucinate a protein binder sequence AlphaFold 
thinks will bind to the target structure. To do this, we minimize PAE and maximize number of contacts at the 
interface and within the binder, and we maximize pLDDT of the binder.
```python
model = mk_afdesign_model(protocol="binder")
model.prep_inputs(pdb_filename="4MZK.pdb", chain="A", binder_len=19)
model.design_3stage(100, 100, 10)
```
Instead of hallucination, you can redesign an existing binder:
```python
model.prep_inputs(pdb_filename="4MZK.pdb", chain="A", binder_chain="T")
```

### partial hallucination
If you have a motif (binding motif, or functional motif) and you want to hallucinate a new scaffold around it,
you can use partial hallucination.
```python
model = mk_afdesign_model(protocol="partial")
model.prep_inputs(pdb_filename="4MZK.pdb", chain="A", pos="1-10,11,20-25")
# TODO
```

# FAQ
#### How do I fixed the FileNotFoundError error?
By default `mk_afdesign_model()` assumes alphafold "params" are saved in the run directory (`data_dir="."`). To override:
```python
model = mk_afdesign_model(..., data_dir="/location/of")
```

#### Can I reuse the same model without needing to recompile?
```python
model.restart()
```
#### How do I change the loss weights?
This can be done using the provided function:
```python
model.set_weights(pae=0.0, plddt=1.0)
```
or the dictionary directly:
```python
model.opt["weights"]["pae"] = 0.0
```
#### How do I control number of recycles used during design?
```python 
model = mk_afdesign_model(num_recycles=1, recycle_mode="average")
# if recycle_mode in ["average","last","sample"] the number of recycles can change during optimization
model.set_opt(recycles=1)
```
- `num_recycles` - number of recycles to use during design (for denovo proteins we find 0 is often enough)
- `recycle_mode` - optimizing across all recycles can be tricky, we experiment with a couple of ways:
  - *last* - use loss from last recycle. (Not recommended, unless you increase number optimization)
  - *sample* - Same as *last* but each iteration a different number of recycles are used. (Previous default).
  - *average* - compute loss at each recycle and average gradients. (Default; Recommended).
  - *add_prev* - average the outputs (dgram, plddt, pae) across all recycles before computing loss.
  - *backprop* - use loss from last recycle, but backprop through all recycles.

#### How do I control which model params are used during design?
By default all five models are used during optimization. If `num_models` > 1, then multiple params are evaluated at each iteration 
and the gradients/losses are averaged. Each iteration a random set of model params are used unless `sample_models=False`.
```python
model = mk_afdesign_model(num_models=1, sample_models=True)
```
- `num_models` - number of model params to use at each iteration.
- `sample_models`:
  - *True* - randomly select models params to use. (Recommended)
  - *False* - use the same model params each iteration.
#### Can I use OpenFold model params for design instead of AlphaFold?
```python
model = mk_afdesign_model(use_openfold=True, use_alphafold=False)
# OR
model.set_opt(use_openfold=True, use_alphafold=False)
```
#### How is contact defined? How do I change it?
By default, 2 [con]tacts per positions are optimized to be within cβ-cβ < 14.0Å and sequence seperation ≥ 9. This can be changed with:
```python
model.set_opt(con=dict(cutoff=8, seqsep=5, num=1))
```
For interface:
```python
model.set_opt(i_con=dict(...))
```
#### For binder hallucination, can I specify the site I want to bind?
```python
model.prep_inputs(..., hotspot="1-10,15,3")
```
#### Can I input more than one chain?
```python
model.prep_inputs(..., chain="A,B")
```
#### Can I design homo-oligomers?
```python
model.prep_inputs(..., copies=2)
# specify interface specific contact and/or pae loss
model.set_weights(i_con=1, i_pae=0)
```
#### For fixed backbone design, how do I force the sequence to be the same for homo-dimer optimization?
```python
model.prep_inputs(pdb_filename="6Q40.pdb", chain="A,B", copies=2, homooligomer=True)
```
WARNING, this functionality assumes the input chains are of equal length.
#### How do I disable certain amino acids?
```python
model.restart(rm_aa="C,W")
```
#### How do I set the random seed for reproducibility?
```python
model.restart(seed=0)
```
#### What are all the different `design_???` methods?
- For **design** we provide 5 different functions:
  - `design_logits()` - optimize *logits* inputs (continious)
  - `design_soft()` - optimize *softmax(logits)* inputs (probabilities)
  - `design_hard()` - optimize *one_hot(logits)* inputs (discrete)

- For complex topologies, we find directly optimizing one_hot encoded sequence `design_hard()` to be very challenging. 
To get around this problem, we propose optimizing in 2 or 3 stages.
  - `design_2stage()` - *soft* → *hard*
  - `design_3stage()` - *logits* → *soft* → *hard*
#### What are all the different losses being optimized?
- general losses
  - *pae*       - minimizes the predicted alignment error
  - *plddt*     - maximizes the predicted LDDT
  - *msa_ent*   - minimize entropy for MSA design (see example at the end of notebook)
  - *pae* and *plddt* values are between 0 and 1 (where lower is better for both)

- fixbb specific losses
  - *dgram_cce* - minimizes the categorical-crossentropy between predicted distogram and one extracted from pdb.
  - *fape*      - minimize difference between coordinates (see AlphaFold paper)
  - we find *dgram_cce* loss to be more stable for design (compared to *fape*)

- hallucination specific losses
  - *con*       - maximize number of contacts. (We find just minimizing *plddt* results in single long helix, 
and maximizing *pae* results in a two helix bundle. To encourage compact structures we add a `con` term)

- binder specific losses
  - *i_pae* - minimize PAE interface of the proteins
  - *pae* - minimize PAE within binder
  - *i_con* - maximize number of contacts at the interface of the proteins
  - *con* - maximize number of contacts within binder

- partial hallucination specific losses
  - *sc_fape* - sidechain-specific fape

# Advanced FAQ
#### loss during Gradient descent is too jumpy, can I do some kind of greedy search towards the end?
Gradient descent updates multiple positions each iteration, which can be a little too aggressive during hard (discrete) mode.
Instead, one can try (`tries`) a few random mutations and accept one with lowest loss. If `use_plddt=True` the random mutations will be biased towards positions with low pLDDT.
```python
model.design_3stage(hard_iters=0)
# set number of model params to evaluate at each iteration
num_models = 2 if model.args["use_templates"] else 5
model.set_opt(models=num_models)
model.design_semigreedy(iters=10, tries=20, use_plddt=True)
```
#### I was getting better results before the major update (19June2022), how do I revert back to the old settings?
We are actively trying to find the best weights `model.opt["weights"]`, settings `model.opt` for each protocol.
Please send us a note if you find something better! To revert back to old settings do this after prepping the model:
- fixbb:
```python
model.restart()
model.set_weights(dgram_cce=1,pae=0.1,plddt=0.1)
model.design_3stage()
```
- hallucination:
```python
model.restart(mode="gumbel")
model.set_weights(pae=1,plddt=1,con=0.5)
model.set_opt(con=dict(binary=True, cutoff=21.6875, num=model._len, seqsep=0))
model.design_2stage(100,100,10)
```
- binder hallucination:
```python
model.restart()
model.set_weights(plddt=0.1, pae=0.1, i_pae=1.0, con=0.1, i_con=0.5)
model.set_opt(con=dict(binary=True, cutoff=21.6875, num=model._binder_len, seqsep=0))
model.set_opt(i_con=dict(binary=True, cutoff=21.6875, num=model._binder_len))
model.design_3stage(100,100,10)
```
#### I don't like your design_??? function, can I write my own with more detailed control?
```python
def design_custom(self):
  # set options
  self.set_opt(dropout=True, soft=0, hard=False)
  # set number of recycles
  self.set_opt(recycles=0)
  # take 100 steps
  for _ in range(100): self._step()
  # increase weight for plddt
  self.set_weights(plddt=2.0)
  # take another 100 steps
  for _ in range(100): self._step()
  # increase number of recycles
  self.set_opt(recycles=1)
  # take another 100 steps
  for _ in range(100): self._step()
  # etc...
  
model = mk_afdesign_model()
design_custom(model)
```
