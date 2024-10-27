import os,sys

from colabdesign.mpnn import mk_mpnn_model
from colabdesign.af import mk_af_model
from colabdesign.shared.protein import pdb_to_string
from colabdesign.shared.parse_args import parse_args

import pandas as pd
import numpy as np
from string import ascii_uppercase, ascii_lowercase
alphabet_list = list(ascii_uppercase+ascii_lowercase)

def get_info(contig):
  F = []
  free_chain = False
  fixed_chain = False
  sub_contigs = [x.split("-") for x in contig.split("/")]
  for n,(a,b) in enumerate(sub_contigs):
    if a[0].isalpha():
      L = int(b)-int(a[1:]) + 1
      F += [1] * L
      fixed_chain = True
    else:
      L = int(b)
      F += [0] * L
      free_chain = True
  return F,[fixed_chain,free_chain]

def main(argv):
  ag = parse_args()
  ag.txt("-------------------------------------------------------------------------------------")
  ag.txt("Designability Test")
  ag.txt("-------------------------------------------------------------------------------------")
  ag.txt("REQUIRED")
  ag.txt("-------------------------------------------------------------------------------------")
  ag.add(["pdb="          ],  None,   str, ["input pdb"])
  ag.add(["loc="          ],  None,   str, ["location to save results"])
  ag.add(["contigs="      ],  None,   str, ["contig definition"])
  ag.txt("-------------------------------------------------------------------------------------")
  ag.txt("OPTIONAL")
  ag.txt("-------------------------------------------------------------------------------------")
  ag.add(["copies="       ],         1,    int, ["number of repeating copies"])
  ag.add(["num_seqs="     ],         8,    int, ["number of mpnn designs to evaluate"])
  ag.add(["initial_guess" ],     False,   None, ["initialize previous coordinates"])
  ag.add(["use_multimer"  ],     False,   None, ["use alphafold_multimer_v3"])
  ag.add(["use_soluble"   ],     False,   None, ["use solubleMPNN"])
  ag.add(["num_recycles=" ],         3,    int, ["number of recycles"])
  ag.add(["rm_aa="],               "C",    str, ["disable specific amino acids from being sampled"])
  ag.add(["num_designs="  ],         1,    int, ["number of designs to evaluate"])
  ag.add(["mpnn_sampling_temp=" ], 0.1,  float, ["sampling temperature used by proteinMPNN"])
  ag.txt("-------------------------------------------------------------------------------------")
  o = ag.parse(argv)

  if None in [o.pdb, o.loc, o.contigs]:
    ag.usage("Missing Required Arguments")

  if o.rm_aa == "":
    o.rm_aa = None

  # filter contig input
  contigs = []
  for contig_str in o.contigs.replace(" ",":").replace(",",":").split(":"):
    if len(contig_str) > 0:
      contig = []
      for x in contig_str.split("/"):
        if x != "0": contig.append(x)
      contigs.append("/".join(contig))

  chains = alphabet_list[:len(contigs)]
  info = [get_info(x) for x in contigs]
  fixed_pos = []
  fixed_chains = []
  free_chains = []
  both_chains = []
  for pos,(fixed_chain,free_chain) in info:
    fixed_pos += pos
    fixed_chains += [fixed_chain and not free_chain]
    free_chains += [free_chain and not fixed_chain]
    both_chains += [fixed_chain and free_chain]

  flags = {"initial_guess":o.initial_guess,
           "best_metric":"rmsd",
           "use_multimer":o.use_multimer,
           "model_names":["model_1_multimer_v3" if o.use_multimer else "model_1_ptm"]}

  if sum(both_chains) == 0 and sum(fixed_chains) > 0 and sum(free_chains) > 0:
    protocol = "binder"
    print("protocol=binder")
    target_chains = []
    binder_chains = []
    for n,x in enumerate(fixed_chains):
      if x: target_chains.append(chains[n])
      else: binder_chains.append(chains[n])
    af_model = mk_af_model(protocol="binder",**flags)
    prep_flags = {"target_chain":",".join(target_chains),
                  "binder_chain":",".join(binder_chains),
                  "rm_aa":o.rm_aa}
    opt_extra = {}
  
  elif sum(fixed_pos) > 0:
    protocol = "partial"
    print("protocol=partial")
    af_model = mk_af_model(protocol="fixbb",
                           use_templates=True,
                           **flags)
    rm_template = np.array(fixed_pos) == 0
    prep_flags = {"chain":",".join(chains),
                  "rm_template":rm_template,
                  "rm_template_seq":rm_template,
                  "copies":o.copies,
                  "homooligomer":o.copies>1,
                  "rm_aa":o.rm_aa}
  else:
    protocol = "fixbb"
    print("protocol=fixbb")
    af_model = mk_af_model(protocol="fixbb",**flags)
    prep_flags = {"chain":",".join(chains),
                  "copies":o.copies,
                  "homooligomer":o.copies>1,
                  "rm_aa":o.rm_aa}

  batch_size = 8
  if o.num_seqs < batch_size:    
    batch_size = o.num_seqs
  
  print("running proteinMPNN...")
  sampling_temp = o.mpnn_sampling_temp
  mpnn_model = mk_mpnn_model(weights="soluble" if o.use_soluble else "original")
  outs = []
  pdbs = []
  for m in range(o.num_designs):
    if o.num_designs == 0:
      pdb_filename = o.pdb
    else:
      pdb_filename = o.pdb.replace("_0.pdb",f"_{m}.pdb")
    pdbs.append(pdb_filename)
    af_model.prep_inputs(pdb_filename, **prep_flags)
    if protocol == "partial":
      p = np.where(fixed_pos)[0]
      af_model.opt["fix_pos"] = p[p < af_model._len]

    mpnn_model.get_af_inputs(af_model)
    outs.append(mpnn_model.sample(num=o.num_seqs//batch_size, batch=batch_size, temperature=sampling_temp))

  if protocol == "binder":
    af_terms = ["plddt","i_ptm","i_pae","rmsd"]
  elif o.copies > 1:
    af_terms = ["plddt","ptm","i_ptm","pae","i_pae","rmsd"]
  else:
    af_terms = ["plddt","ptm","pae","rmsd"]

  labels = ["design","n","score"] + af_terms + ["seq"]
  data = []
  best = {"rmsd":np.inf,"design":0,"n":0}
  print("running AlphaFold...")
  os.system(f"mkdir -p {o.loc}/all_pdb")
  with open(f"{o.loc}/design.fasta","w") as fasta:
    for m,(out,pdb_filename) in enumerate(zip(outs,pdbs)):
      out["design"] = []
      out["n"] = []
      af_model.prep_inputs(pdb_filename, **prep_flags)
      for k in af_terms: out[k] = []
      for n in range(o.num_seqs):
        out["design"].append(m)
        out["n"].append(n)
        sub_seq = out["seq"][n].replace("/","")[-af_model._len:]
        af_model.predict(seq=sub_seq, num_recycles=o.num_recycles, verbose=False)
        for t in af_terms: out[t].append(af_model.aux["log"][t])
        if "i_pae" in out:
          out["i_pae"][-1] = out["i_pae"][-1] * 31
        if "pae" in out:
          out["pae"][-1] = out["pae"][-1] * 31
        rmsd = out["rmsd"][-1]
        if rmsd < best["rmsd"]:
          best = {"design":m,"n":n,"rmsd":rmsd}
        af_model.save_current_pdb(f"{o.loc}/all_pdb/design{m}_n{n}.pdb")
        af_model._save_results(save_best=True, verbose=False)
        af_model._k += 1
        score_line = [f'design:{m} n:{n}',f'mpnn:{out["score"][n]:.3f}']
        for t in af_terms:
          score_line.append(f'{t}:{out[t][n]:.3f}')
        print(" ".join(score_line)+" "+out["seq"][n])
        line = f'>{"|".join(score_line)}\n{out["seq"][n]}'
        fasta.write(line+"\n")
      data += [[out[k][n] for k in labels] for n in range(o.num_seqs)]
      af_model.save_pdb(f"{o.loc}/best_design{m}.pdb")

  # save best
  with open(f"{o.loc}/best.pdb", "w") as handle:
    remark_text = f"design {best['design']} N {best['n']} RMSD {best['rmsd']:.3f}"
    handle.write(f"REMARK 001 {remark_text}\n")
    handle.write(open(f"{o.loc}/best_design{best['design']}.pdb", "r").read())
    
  labels[2] = "mpnn"
  df = pd.DataFrame(data, columns=labels)
  df.to_csv(f'{o.loc}/mpnn_results.csv')

if __name__ == "__main__":
   main(sys.argv[1:])
