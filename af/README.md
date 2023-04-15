# AfDesign (v1.1.1)
### Google Colab
<a href="https://colab.research.google.com/github/sokrypton/ColabDesign/blob/v1.1.1/af/design.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

# Updates
- Jump to [Previous Updates](#previous-updates)
- **15Oct2022** - v1.1.0
  - integrating proteinMPNN!
  - bugfix for sidechain loss
- **17Nov2022**
  - updating pae/plddt loss calculation to be consistent with pae/plddt outputs
- **24Dec2022** - v1.1.1
  - adding af_pseudo_diffusion examples
  - updating to alphafold-multimer v2.3.0
  - enabling fused_triangle_multiplication by default
- **21Jan2023**
  - add support for bfloat16 (enabled by default
- **01Mar2023**
  - adding support for [RfDiffusion](/rf)
### setup
first install jax (with GPU support)
```bash
pip install "jax[cuda]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```
second install colabdesign
```bash
pip install git+https://github.com/sokrypton/ColabDesign.git@v1.1.1

# download alphafold weights
mkdir params
curl -fsSL https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar | tar x -C params
```
By default `mk_afdesign_model()` assumes alphafold "params" are saved in the run directory (`data_dir="."`). To override:
```python
model = mk_afdesign_model(..., data_dir="/location/of")
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
model.set_seq(mode="gumbel")
model.design_soft(50)
model.set_seq(model.aux["seq"]["pseudo"])
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
you can use partial hallucination. Or you have a protein and you want to extend one of the loops.
```python
af_model = mk_afdesign_model(protocol="partial")
af_model.prep_inputs(pdb_filename="6MRR.pdb", chain="A", pos="3-30,33-68", length=100)
af_model.rewire(loops=[36])
```
# FAQ

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
# if recycle_mode in ["average",last","sample","first"] the number of recycles can change during optimization
model.set_opt(num_recycles=1)
```
- `num_recycles` - number of recycles to use during design (for denovo proteins we find 0 is often enough)
- `recycle_mode` - optimizing across all recycles can be tricky, we experiment with a couple of ways:
  - *last* - use loss from last recycle. (Default)
  - *average* - compute loss at each recycle and average gradients. (Previous default from v.1.0.5)
  - *sample* - Same as *last* but each iteration a different number of recycles are used.
  - *first* - use loss from first recycle.
  - *add_prev* - average the outputs (dgram, plddt, pae) across all recycles before computing loss.
  - *backprop* - use loss from last recycle, but backprop through all recycles.

#### How do I control which model params are used during design?
By default all five models are used during optimization. If `num_models` > 1, then multiple params are evaluated at each iteration and the gradients/losses are averaged. Each iteration a random set of model params are used unless `sample_models=False`.
```python
model = mk_afdesign_model(num_models=1, sample_models=True)
# or
model.set_opt(num_models=1, sample_models=True)
```
- `num_models` - number of model params to use at each iteration.
- `sample_models`:
  - *True* - randomly select models params to use. (Recommended)
  - *False* - use the same model params each iteration.
You can also specify exactly which models are used during any of the design protocols:
```python
model.design_(num_models=1, sample_models=True, models=[0,2,3])
# or
model.design_(num_models=2, sample_models=False, models=["model_1_ptm","model_3_ptm"])
```
#### Can I use OpenFold model params for design instead of AlphaFold?
You may need to download them:
```bash
  for W in openfold_model_ptm_1 openfold_model_ptm_2 openfold_model_no_templ_ptm_1
  do wget -qnc https://files.ipd.uw.edu/krypton/openfold/${W}.npz -P params; done
```
Once downloaded:
```python
model = mk_afdesign_model(use_openfold=True, use_alphafold=False)
```
#### For binder hallucination, can I specify the site I want to bind?
```python
model.prep_inputs(..., hotspot="1-10,15,3")
```
#### Can I input more than one chain?
```python
model.prep_inputs(..., chain="A,B")
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
To get around this problem, we propose optimizing in 3 stages or first learning logits then switching to semigreedy optimization.
  - `design_3stage()` - gradient based optimization (GD) (logits → soft → hard)
  - `design_semigreedy(tries=X)` - tries X random mutations, accepts those that decrease loss
  - `design_pssm_semigreey(tries=X)` - uses GD to get a sequence profile (PSSM), then uses the PSSM to bias semigreedy opt. (Recommended)

#### What are all the different losses being optimized?
- general losses
  - *pae*       - minimizes the predicted alignment error
  - *plddt*     - maximizes the predicted LDDT
  - *pae* and *plddt* values are between 0 and 1 (where lower is better for both)

- fixbb specific losses
  - *dgram_cce* - minimizes the categorical-crossentropy between predicted distogram and one extracted from pdb.
  - *fape*      - minimize difference between coordinates (see AlphaFold paper)
  - we find *dgram_cce* loss to be more stable for design (compared to *fape*)

- hallucination specific losses
  - *con*       - maximize `1` contacts per position. `model.set_opt("con",num=1)`

- binder specific losses
  - *pae* - minimize PAE at interface and within binder
  - *con* - - maximize `2` contacts per binder position, within binder. `model.set_opt("con",num=2)`
  - *i_con* - maximize `1` contacts per binder position `model.set_opt("i_con",num=1)`

- partial hallucination specific losses
  - *sc_fape* - sidechain-specific fape

#### How is contact defined? How do I change it?
By default, 2 [con]tacts per positions are optimized to be within cβ-cβ < 14.0Å and sequence seperation ≥ 9. This can be changed with:
```python
model.set_opt(con=dict(cutoff=8, seqsep=5, num=1))
```
For interface:
```python
model.set_opt(i_con=dict(...))
```

#### Optax Optimizers
By default, we use stochastic gradient descent `set_optimizer(optimizer="sgd", learning_rate=0.1, norm_seq_grad=True)` for optimization. This seems to work quite well for the default problems. But if you want to try other optimizers, ColabDesign is now fully integrated with all [Optax optimizers](https://optax.readthedocs.io/en/latest/api.html).
Example how to change optimizer. Note the default learning_rate of 0.1 was calibrated for sgd with gradient normalization. A different learning_rate may be more optimal for other optimization settings.
```python
model = mk_afdesign_model(optimizer="adam", learning_rate=0.01)
```
Or for more control (or to change settings after model initialization), use:
```python
model.set_optimizer(optimizer="adam", learning_rate=0.01, b1=0.9, b2=0.999)
```
By default, the gradients for the sequence parameters are normalized. We find this helps with convergence. To disable:
```python
model.set_opt(norm_seq_grad=False)
```

# Advanced FAQ
#### loss during Gradient descent is too jumpy, can I do some kind of greedy search towards the end?
Gradient descent updates multiple positions each iteration, which can be a little too aggressive during hard (discrete) mode.
Instead, one can try (`tries`) a few random mutations and accept one with lowest loss. If `use_plddt=True` the random mutations will be biased towards positions with low pLDDT.
```python
model.design_3stage(hard_iters=0)
# set number of model params to evaluate at each iteration
num_models = 2 if model.args["use_templates"] else 5
model.design_semigreedy(iters=10, tries=20, num_models=num_models, use_plddt=True)
```
#### I was getting better results before the major update (19June2022), how do I revert back to the old settings?
We are actively trying to find the best weights `model.opt["weights"]`, settings `model.opt` for each protocol.
Please send us a note if you find something better! To revert back to old settings do this after prepping the model:
- fixbb:
```python
model.set_weights(dgram_cce=1, pae=0.1, plddt=0.1)
model.design_3stage()
```
- hallucination:
```python
model.set_seq(mode="gumbel")
model.set_weights(pae=1, plddt=1, con=0.5)
model.set_opt("con", binary=True, cutoff=21.6875, num=model._len, seqsep=0)
model.design_2stage(100, 100, 10)
```
- binder hallucination:
```python
model.set_weights(plddt=0.1, pae=0.1, i_pae=1.0, con=0.1, i_con=0.5)
model.set_opt("con", binary=True, cutoff=21.6875, num=model._binder_len, seqsep=0)
model.set_opt("i_con", binary=True, cutoff=21.6875, num=model._target_len)
model.design_3stage(100, 100, 10)
```
#### I don't like your design_??? function, can I write my own with more detailed control?
```python
def design_custom(self):
  # set options
  self.set_opt(dropout=True, soft=False, hard=False)
  # set number of recycles
  self.set_opt(num_recycles=0)
  # take 100 steps
  for _ in range(100): self.step()
  # increase weight for plddt
  self.set_weights(plddt=2.0)
  # take another 100 steps
  for _ in range(100): self.step()
  # increase number of recycles
  self.set_opt(num_recycles=1)
  # take another 100 steps
  for _ in range(100): self.step()
  # etc...
  
model = mk_afdesign_model()
design_custom(model)
```

#### custom callback examples
Looking for more control over afdesign? The callback functions have gotten much smarter. Based on your input arguments, it will automatically fetch the variable of interest. You can now define your own custom losses, params to optimize, modify inputs before alphafold is run, and modify auxiliary outputs.
```python
def custom_pre_callback(inputs, aux, opt, key):
  inputs["aatype"] = inputs["aatype"].at[:].set(0)
  aux["pre"] = opt["pre"] + jax.random.randint(key,[],0,10)

def custom_post_callback(outputs, aux):
  aux["post"] = outputs["structure_module"]

def custom_loss_callback(outputs, params):
  loss = jnp.square(outputs["structure_module"]["final_atom14_positions"] + params["custom_param"]).mean()
  return {"custom_loss":loss}

af_model = mk_afdesign_model(protocol="fixbb",
                             pre_callback=custom_pre_callback,
                             post_callback=custom_post_callback,
                             loss_callback=custom_loss_callback)
af_model._params["custom_param"] = 1.0
af_model.opt["weights"]["custom_loss"] = 0.1
af_model.opt["pre"] = 100

af_model.prep_inputs(pdb_filename=get_pdb("1TEN"), chain="A")
```

# Previous Updates
- **24Feb2022** - "Beta" branch started. Refactoring code to allow homooligomeric hallucination/design and averaging gradients across recycles (which is now the default).
Minor changes changes include renaming intra_pae/inter_con to pae/con and inter_pae/inter_con to i_pae/i_con for clarity.
- **28Feb2022** - We find backprop through structure module to be unstable, all functions have been updated to only use distogram by default. The definition of contact has changed to minimize entropy within distance cutoff.
- **02May2022** - The `design.py` code has been split up into multiple python files under `src/`
- **14May2022** - Adding support for partial hallucination (if you want to constrain one part and generate structure/sequence for rest).
- **19June2022** - "Beta" branch is now the "Main" branch. WARNING: Lots of default settings and weights were changed. [Click here](#i-was-getting-better-results-before-the-major-update-19june2022-how-do-i-revert-back-to-the-old-settings) for info on how to revert back to old settings. 
- **28June2022** - v1.0.1 - Major code reorganization/refactoring to add support for callbacks (to allow integration w/ other tools during design) and to avoid clashes with existing trrosetta/alphafold installations. (eg. `af → colabdesign`, `af.src → colabdesign.af` and `alphafold → colabdesign.af.alphafold`).
- **05July2022** - v1.0.2 - Major code cleanup, removing duplicate code. Adding support for custom loss functions.
- **11July2022** - v1.0.3 - Improved homo-oligomeric support. RMSD and dgram losses have been refactored to automatically save aligned coordinates. Multimeric coordinates now saved with chain identifiers.
- **23July2022** - v1.0.4 - Adding support for openfold weights. To enable set `mk_afdesign_model(..., use_openfold=True)`.
- **31July2022** - v1.0.5 - Refactoring to add support for swapping batch features without recompile. Allowing for implementation of [AF2Rank](https://github.com/sokrypton/ColabDesign/blob/main/af/examples/AF2Rank.ipynb)!
- **09Sept2022** - v1.0.6
  - support for alphafold-multimer `model = mk_afdesign_model(..., use_multimer=True)`
  - support for experimentally resolved loss `model.set_weights(exp_res=1)`
  - support for multichain design/hallucination for fixbb, hallucination and partial protocols: `model.prep_inputs(..., copies=2)`
  - support to fix the sequence for certain positions `model.prep_inputs(..., fix_pos="1-10")` (supported in protocols "fixbb" and "partial")
  - binder protocol improved, prior protocol would try to optimize number of contacts per target, new default is to optimize number of contacts per binder position. Number of contacts per binder position can be controlled with `model.set_opt("i_con",num=1)` and number of positions that should be contact with `model.set_opt("i_con",num_pos=5)`
  - implementing David Jones'-like protocol for semi-greedy optimization, where positions are selected based on plddt, and after 20 tries, the mutation that decreasing loss the most is accepted. `model.design_semigreedy()`
  - WARNING: the returned pLDDT is now in the "correct" direction (higher is better)
  - removing recycle dimension from the input features (to standardize with multimer inputs)
  - removing all dependence on TensorFlow
- **14Sept2022** - v1.0.7
  - refactoring design.py to add `design_pssm_semigreedy()` protocol, which is a wrapper around `design_semigreedy(seq_logits=)`, and can be used to input/learn PSSM for biased optimization.
  - adding example [peptide_binder_design.ipynb](https://colab.research.google.com/github/sokrypton/ColabDesign/blob/main/af/examples/peptide_binder_design.ipynb) targeted for peptide binder hallucination/design.
  - adding [finer control](#how-do-i-control-which-model-params-are-used-during-design) over what models are used during optimization.
  - fixing RAM memory leaks, `clear_mem()` now also does garbage collection
  - fixing integration with TrDesign that got broken in v1.0.6
- **22Sept2022** - v1.0.8
  - [custom callback functions](#custom-callback-examples) (\[pre|loss|pos\]_callback) have been refactored to be more flexible.
    - Supported input arguments include: ["inputs", "outputs", "params", "opt", "seq", "aux", "key"]. 
    - The pre_callback function can be used to modify inputs before prediction, loss_callback to add cutstom loss.
  - adding support for [Optax optimizers](#optax-optimizers)
- **24Sept2022** - v1.0.9
  - adding [contrib section](/af/contrib) where user contributed modifications and protocols will go.
