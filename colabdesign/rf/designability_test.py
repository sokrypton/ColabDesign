import os,sys
from colabdesign.mpnn import mk_mpnn_model
from colabdesign.af import mk_af_model
from colabdesign.shared.protein import pdb_to_string
import pandas as pd
import numpy as np
from string import ascii_uppercase, ascii_lowercase
alphabet_list = list(ascii_uppercase+ascii_lowercase)

_,pdb,loc,contigs,copies,initial_guess = sys.argv
copies = int(copies)
initial_guess = bool(initial_guess)

def get_info(contig):
  F = []
  fixed_chain = True
  sub_contigs = [x.split("-") for x in contig.split("/")]
  for n,(a,b) in enumerate(sub_contigs):
    if a[0].isalpha():
      L = int(b)-int(a[1:]) + 1
      F += [1] * L
    else:
      L = int(b)
      F += [0] * L
      fixed_chain = False
  return F,fixed_chain

contigs = contigs.split(":")
chains = alphabet_list[:len(contigs)]
info = [get_info(x) for x in contigs]
fixed_chains = [y for x,y in info]
fixed_pos = sum([x for x,y in info],[])

if sum(fixed_chains) > 0 and sum(fixed_chains) < len(fixed_chains):
  protocol = "binder"
  print("protocol=binder")
  target_chains = []
  binder_chains = []
  for n,x in enumerate(fixed_chains):
    if x: target_chains.append(chains[n])
    else: binder_chains.append(chains[n])
  af_model = mk_af_model(protocol="binder",
                         model_names=["model_1_ptm"],
                         best_metric="rmsd",
                         initial_guess=initial_guess)
  af_model.prep_inputs(pdb,
                       target_chain=",".join(target_chains),
                       binder_chain=",".join(binder_chains))
elif sum(fixed_pos) > 0:
  protocol = "partial"
  print("protocol=partial")
  af_model = mk_af_model(protocol="fixbb",
                         model_names=["model_1_ptm"],
                         use_templates=True,
                         best_metric="rmsd",
                         initial_guess=initial_guess)
  rm_template = np.array(fixed_pos) == 0
  af_model.prep_inputs(pdb,
                       chain=",".join(chains),
                       rm_template=rm_template,
                       rm_template_seq=rm_template,
                       copies=copies,
                       homooligomer=copies>1)
  p = np.where(fixed_pos)[0]
  af_model.opt["fix_pos"] = p[p < af_model._len]

else:
  protocol = "fixbb"
  print("protocol=fixbb")
  af_model = mk_af_model(protocol="fixbb",
                         model_names=["model_4_ptm"], 
                         best_metric="rmsd",
                         initial_guess=initial_guess)
  af_model.prep_inputs(pdb,
                       chain=",".join(chains),
                       copies=copies,
                       homooligomer=copies>1)

print("running proteinMPNN...")
num_seqs = 8
sampling_temp = 0.1
mpnn_model = mk_mpnn_model()
mpnn_model.get_af_inputs(af_model)
out = mpnn_model.sample(num=num_seqs//8, batch=8, temperature=sampling_temp)
print("running AlphaFold...")
if protocol == "binder":
  af_terms = ["plddt","i_ptm","i_pae","rmsd"]
elif copies > 1:
  af_terms = ["plddt","ptm","i_ptm","pae","i_pae","rmsd"]
else:
  af_terms = ["plddt","ptm","pae","rmsd"]
for k in af_terms: out[k] = []
os.system(f"mkdir -p {loc}/all_pdb")
with open(f"{loc}/design.fasta","w") as fasta:
  for n in range(num_seqs):
    seq = out["seq"][n][-af_model._len:]
    af_model.predict(seq=seq, num_recycles=3, verbose=False)
    for t in af_terms: out[t].append(af_model.aux["log"][t])
    if "i_pae" in out:
      out["i_pae"][-1] = out["i_pae"][-1] * 31
    if "pae" in out:
      out["pae"][-1] = out["pae"][-1] * 31
      
    af_model.save_current_pdb(f"{loc}/all_pdb/n{n}.pdb")
    af_model._save_results(save_best=True, verbose=False)
    af_model._k += 1
    score_line = [f'mpnn:{out["score"][n]:.3f}']
    for t in af_terms:
      score_line.append(f'{t}:{out[t][n]:.3f}')
    print(n, " ".join(score_line)+" "+seq)
    line = f'>{"|".join(score_line)}\n{seq}'
    fasta.write(line+"\n")

af_model.save_pdb(f"{loc}/best.pdb")

labels = ["score"] + af_terms + ["seq"]
data = [[out[k][n] for k in labels] for n in range(num_seqs)]
labels[0] = "mpnn"

df = pd.DataFrame(data, columns=labels)
df.to_csv(f'{loc}/mpnn_results.csv')