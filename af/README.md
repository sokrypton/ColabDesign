# AfDesign (beta branch)
### Google Colab
<a href="https://colab.research.google.com/github/sokrypton/ColabDesign/blob/beta/af/design.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

### setup
```bash
git clone https://github.com/sokrypton/af_backprop.git
pip -q install biopython dm-haiku==0.0.5 ml-collections py3Dmol
mkdir params
curl -fsSL https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar | tar x -C params
wget -qnc https://raw.githubusercontent.com/sokrypton/ColabFold/main/beta/colabfold.py
wget -qnc https://raw.githubusercontent.com/sokrypton/ColabDesign/beta/af/design.py
```
```python
import numpy as np
from IPython.display import HTML
from design import mk_design_model, clear_mem
import sys
sys.path.append('af_backprop')
```
### fixed backbone design
For a given protein backbone, generate/design a new sequence that AlphaFold thinks folds into that conformation
```python
model = mk_design_model(protocol="fixbb")
model.prep_inputs(pdb_filename="1TEN.pdb", chain="A")
model.design_3stage()
```
### hallucination
For a given length, generate/hallucinate a protein sequence that AlphaFold thinks folds into a well structured 
protein (high plddt, low pae, many contacts).
```python
model = mk_design_model(protocol="hallucination")
model.prep_inputs(length=100, seq_init="gumbel")
model.design_2stage()
```
### binder hallucination
For a given protein target and protein binder length, generate/hallucinate a protein binder sequence AlphaFold 
thinks will bind to the target structure. To do this, we minimize PAE and maximize number of contacts at the 
interface and within the binder, and we maximize pLDDT of the binder.
```python
model = mk_design_model(protocol="binder")
model.prep_inputs(pdb_filename="4MZK.pdb", chain="A", binder_len=19)
model.design_3stage(soft_iters=100, temp_iters=100, hard_iters=10)
```
# Updates
- **24Feb2022** - Refactoring code to allow homooligomeric hallucination/design and averaging gradients across recycles (which is now the default).
Minor changes changes include renaming intra_pae/inter_con to pae/con and inter_pae/inter_con to i_pae/i_con for clarity.
- **28Feb2022** - We find backprop through structure module to be unstable, all functions have been updated to only use distogram by default.
# FAQ
#### Can I reuse the same model without needing to recompile?
```python
model.restart()
```
#### How do I change the loss weights?
```python
model.opt["weights"].update({"pae":0.0,"plddt":1.0})
model.opt["weights"]["pae"] = 0.0
```
WARNING: When setting weights be careful to use floats (instead of `1`, use `1.0`), otherwise this triggers recompile.
#### How do I control number of recycles used during design?
```python 
model = mk_design_model(num_recycles=1, recycle_mode="sample")
```
- `num_recycles` - max number of recycles to use during design (for denovo proteins we find 0 is often enough)
- `recycle_model` - When `num_recycles` > 0:
  - *average* - compute loss at each recycle and average gradients. (Recommend).
  - *sample* - at each iteration, randomly select number of recycles to use, use gradients from last recycle.
  - *add_prev* - average the logits (dgram, plddt, pae) across all recycles.
  - *last* - only use gradients from last recycle. (NOT recommended, unless you want to start with 0 and increase later).
#### How do I control which model params are used during design?
By default all five models are used during optimization. If `num_models` > 1, then multiple params are evaluated at each iteration 
and the gradients/losses are averaged. Each iteration a random set of model params are used unless `model_mode="fixed"`.
```python
model = mk_design_model(num_models=1, model_mode="sample", model_parallel=False)
```
- `num_models` - number of model params to use at each iteration.
- `model_mode`:
  - *sample* - randomly select models params to use. (Recommended)
  - *fixed* - use the same model params each iteration.
- `model_parallel` - run model params in parallel if `num_models` > 1. By default, the model params are evaluated in serial,
if you have access to high-end GPU, you can run all model params in parallel by enabling this flag. 
#### How is contact defined? How do I change it?
[con]tact is defined as cβ-cβ < 14.0Å (where distogram bins < 14Å are summed) and sequence seperation ≥ 9. This can be changed with:
```python
model.opt["con"].update({"cutoff":8.0,"seqsep":5})
```
For interface:
```python
model.opt["i_con"].update(...)
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
model.opt["weights"].update({"i_con":0.0, "i_pae":0.0})
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

# Advanced FAQ
#### I don't like your design_??? function, can I write my own with more detailed control?
```python
def design_custom(self):
  # set options
  self.opt.update({"dropout":True, "soft":0.0, "hard":False"})
  # set number of recycles
  self.opt["recycles"] = 0
  # take 100 steps
  for _ in range(100): self._step()
  # increase weight for plddt
  self.opt["weights"].update({"plddt":2.0})
  # take another 100 steps
  for _ in range(100): self._step()
  # increase number of recycles
  self.opt["recycles"] = 1
  # take another 100 steps
  for _ in range(100): self._step()
  # etc...
  
model = mk_design_model()
design_custom(model)
```
