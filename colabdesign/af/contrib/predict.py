import os, time
import re, tempfile
from IPython.display import HTML
from colabdesign.af.inputs import copy_dict
from colabdesign.af.prep import prep_pdb
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.shared.protein import _np_get_cb
from colabdesign.shared.parsers import parse_a3m

import numpy as np
import matplotlib.pyplot as plt
from string import ascii_uppercase, ascii_lowercase
import hashlib
import random
import jax
from collections import OrderedDict, Counter

import Bio
from Bio.Align import substitution_matrices

aligner = Bio.Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

pymol_color_list = ["#33ff33","#00ffff","#ff33cc","#ffff00","#ff9999","#e5e5e5","#7f7fff","#ff7f00",
                    "#7fff7f","#199999","#ff007f","#ffdd5e","#8c3f99","#b2b2b2","#007fff","#c4b200",
                    "#8cb266","#00bfbf","#b27f7f","#fcd1a5","#ff7f7f","#ffbfdd","#7fffff","#ffff7f",
                    "#00ff7f","#337fcc","#d8337f","#bfff3f","#ff7fff","#d8d8ff","#3fffbf","#b78c4c",
                    "#339933","#66b2b2","#ba8c84","#84bf00","#b24c66","#7f7f7f","#3f3fa5","#a5512b"]

alphabet_list = list(ascii_uppercase+ascii_lowercase)
aa_order = {a:n for n,a in enumerate(residue_constants.restypes + ["X","-"])}
order_aa = {n:a for a,n in aa_order.items()}

def get_hash(x):
  return hashlib.sha1(x.encode()).hexdigest()

def get_unique_sequences(seq_list):
  unique_seqs = list(OrderedDict.fromkeys(seq_list))
  return unique_seqs

def parse_hhalign_output(filename):
  with open(filename, "r") as f:
    lines = f.readlines()
  start_indices, sequences = [None,None], ["",""]
  for line in lines:
    parts = line.split()
    if line.startswith('Q '):
      if start_indices[0] is None:
        start_indices[0] = int(parts[2]) - 1
      sequences[0] += parts[3]
    if line.startswith('T '):
      if start_indices[1] is None:
        start_indices[1] = int(parts[2]) - 1
      sequences[1] += parts[3]
  return sequences, start_indices

def get_dgram(positions=None, seq=None,
              dist=None, num_bins=39, min_bin=3.25, max_bin=50.75):
  if dist is None:
    if seq is None:
      atoms = {k:positions[...,residue_constants.atom_order[k],:] for k in ["N","CA","C"]}
      c = _np_get_cb(**atoms, use_jax=False)
    else:
      ca = positions[...,residue_constants.atom_order["CA"],:]
      cb = positions[...,residue_constants.atom_order["CB"],:]
      is_gly = seq==residue_constants.restype_order["G"]
      c = np.where(is_gly[:,None],ca,cb)
    dist = np.sqrt(np.square(c[None,:] - c[:,None]).sum(-1,keepdims=True))
  lower_breaks = np.linspace(min_bin, max_bin, num_bins)
  lower_breaks = lower_breaks
  upper_breaks = np.append(lower_breaks[1:],1e8)
  def get_bins(d):
    return ((d > lower_breaks) * (d < upper_breaks)).astype(float)
  return get_bins(dist)

def get_template_idx(query_sequence, target_sequence, query_a3m=None, target_a3m=None, align_fn=None):

  if align_fn is None:
    align_fn = lambda a,b,**c: list(aligner.align(a[0],b[0])[0]),(0,0)

  if len(query_sequence) == len(target_sequence):
    i = np.arange(len(query_sequence))
    return i,i
  else:
    X,(i,j) = align_fn(query_sequence, target_sequence, query_a3m=query_a3m, target_a3m=target_a3m)
    print(f"Q\t{X[0]}\nT\t{X[1]}")
    idx_i, idx_j = [], []
    for q,t in zip(*X):
      if q != "-" and t != "-":
          idx_i.append(i)
          idx_j.append(j)
      if q != "-": i += 1
      if t != "-": j += 1
    return np.array(idx_i), np.array(idx_j)

def get_template_feats(pdbs, chains, query_seq, query_a3m=None, 
  copies=1, propagate_to_copies=True, use_seq=False, use_dgram=True, 
  get_pdb_fn=None, align_fn=None):

  if get_pdb_fn is None:
    get_pdb_fn = lambda x:x
  def cat(x):
    if len(x) > 1:
      return jax.tree_map(lambda *y:np.concatenate(y,0), *x)
    else:
      return x[0]
  if isinstance(pdbs,str): pdbs = pdbs.split(":")
  if isinstance(chains,str): chains = chains.split(":")
  # load and concatenate pdb(s) and chain(s)
  N,X,L,P,C = [],[],[],[],[]
  for n,(pdb,chain) in enumerate(zip(pdbs,chains)):
    pdb_filename = get_pdb_fn(pdb)
    if isinstance(chain,str): chain = chain.split(",")
    for c in chain:
      info = prep_pdb(pdb_filename, c, ignore_missing=True)
      N.append(n)
      X.append(info)
      L.append(info.pop("lengths")[0])
      P.append(pdb)
      C.append(c)
  X = cat(X)
  
  if use_dgram:
    mask = np.ones((sum(L),sum(L)))
    bL = np.cumsum(L)
    aL = bL - np.array(L)
    for i in range(len(L)):
      for j in range(i+1,len(L)):
        if N[i] != N[j] or (P[i] == P[j] and C[i] == C[j]) or P[i] != P[j]:
          mask[aL[i]:bL[i],aL[j]:bL[j]] = 0
          mask[aL[j]:bL[j],aL[i]:bL[i]] = 0
  
  seq = X["batch"]["aatype"] if use_seq else None
  X["batch"]["dgram"] = get_dgram(X["batch"]["all_atom_positions"], seq) * mask[...,None]
  template_seq = "".join(order_aa.get(a,"X") for a in X["batch"]["aatype"])
  
  Q = query_seq if propagate_to_copies else query_seq * copies
  i, j = get_template_idx(Q, template_seq, query_a3m=query_a3m, align_fn=align_fn)
  
  L = len(Q)
  batch = {"aatype":np.zeros(L,int),
           "all_atom_mask":np.zeros((L,37)),
           "all_atom_positions":np.zeros((L,37,3))}
  if use_dgram: batch["dgram"] = np.zeros((L,L,39))
  
  for k in batch:
    if k == "dgram":
      batch[k][i[:,None,None],i[None,:,None],:] = X["batch"][k][j[:,None,None],j[None,:,None]]
    else:
      batch[k][i] = X["batch"][k][j]

  if propagate_to_copies and copies > 1:
    tile_shape = {
      "aatype":[copies], 
      "all_atom_mask":[copies,1], 
      "all_atom_positions":[copies,1,1]}
    if use_dgram: tile_shape["dgram"] = [copies,copies,1]
    
    for k,s in tile_shape.items():
      batch[k] = np.tile(batch[k],s)

    if use_dgram:
      # mask non-diagonal
      x = batch["dgram"]
      mask = np.eye(copies)[:,None,:,None,None]
      batch["dgram"] = (x.reshape(copies,L,copies,L,-1) * mask).reshape(x.shape)

  return batch

def plot_msa(msa, Ls, sort_lines=True, dpi=100):
  seq = msa[0]
  Ln = np.cumsum([0] + Ls)
  gap = msa != 21
  qid = msa == seq
  gapid = np.stack([gap[:,Ln[i]:Ln[i+1]].max(-1) for i in range(len(Ls))],-1)
  lines = []
  Nn = []
  for g in np.unique(gapid, axis=0):
    i = np.where((gapid == g).all(axis=-1))
    qid_ = qid[i]
    gap_ = gap[i]
    seqid = np.stack([qid_[:,Ln[i]:Ln[i+1]].mean(-1) for i in range(len(Ls))],-1).sum(-1) / (g.sum(-1) + 1e-8)
    non_gaps = gap_.astype(float)
    non_gaps[non_gaps == 0] = np.nan
    if sort_lines:
      lines_ = non_gaps[seqid.argsort()] * seqid[seqid.argsort(),None]
    else:
      lines_ = non_gaps[::-1] * seqid[::-1,None]
    Nn.append(len(lines_))
    lines.append(lines_)
  
  Nn = np.cumsum(np.append(0,Nn))
  lines = np.concatenate(lines,0)
  plt.figure(figsize=(8,5), dpi=dpi)
  plt.title("Sequence coverage")
  plt.imshow(lines,
        interpolation='nearest', aspect='auto',
        cmap="rainbow_r", vmin=0, vmax=1, origin='lower',
        extent=(0, lines.shape[1], 0, lines.shape[0]))
  for i in Ln[1:-1]:
    plt.plot([i,i],[0,lines.shape[0]],color="black")
  for j in Nn[1:-1]:
    plt.plot([0,lines.shape[1]],[j,j],color="black")
  
  plt.plot((np.isnan(lines) == False).sum(0), color='black')
  plt.xlim(0,lines.shape[1])
  plt.ylim(0,lines.shape[0])
  plt.colorbar(label="Sequence identity to query")
  plt.xlabel("Positions")
  plt.ylabel("Sequences")
  return plt

def plot_plddt_legend(dpi=100):
  thresh = ['plDDT:','Very low (<50)','Low (60)','OK (70)','Confident (80)','Very high (>90)']
  plt.figure(figsize=(1,0.1),dpi=dpi)
  ########################################
  for c in ["#FFFFFF","#FF0000","#FFFF00","#00FF00","#00FFFF","#0000FF"]:
    plt.bar(0, 0, color=c)
  plt.legend(thresh, frameon=False,
             loc='center', ncol=6,
             handletextpad=1,
             columnspacing=1,
             markerscale=0.5,)
  plt.axis(False)
  return plt

def plot_ticks(Ls, axes=None):
  if axes is None: axes = plt.gca()
  Ln = sum(Ls)
  L_prev = 0
  for L_i in Ls[:-1]:
    L = L_prev + L_i
    L_prev += L_i
    plt.plot([0,Ln],[L,L],color="black")
    plt.plot([L,L],[0,Ln],color="black")
  ticks = np.cumsum([0]+Ls)
  ticks = (ticks[1:] + ticks[:-1])/2
  axes.set_yticks(ticks)
  axes.set_yticklabels(alphabet_list[:len(ticks)])

def plot_confidence(plddt, pae=None, Ls=None, dpi=100):
  use_ptm = False if pae is None else True
  if use_ptm:
    plt.figure(figsize=(10,3), dpi=dpi)
    plt.subplot(1,2,1);
  else:
    plt.figure(figsize=(5,3), dpi=dpi)
  plt.title('Predicted lDDT')
  plt.plot(plddt)
  if Ls is not None:
    L_prev = 0
    for L_i in Ls[:-1]:
      L = L_prev + L_i
      L_prev += L_i
      plt.plot([L,L],[0,100],color="black")
  plt.ylim(0,100)
  plt.ylabel('plDDT')
  plt.xlabel('position')
  if use_ptm:
    plt.subplot(1,2,2);plt.title('Predicted Aligned Error')
    Ln = pae.shape[0]
    plt.imshow(pae,cmap="bwr",vmin=0,vmax=30,extent=(0, Ln, Ln, 0))
    if Ls is not None and len(Ls) > 1: plot_ticks(Ls)
    plt.colorbar()
    plt.xlabel('Scored residue')
    plt.ylabel('Aligned residue')
  return plt

def get_msa(seqs, jobname, cov=50, id=90, qid=10, max_msa=4096, mode="unpaired_paired", 
  mmseqs2_fn=None, hhfilter_fn=None, verbose=True, output_a3m=None,
  do_not_filter=False, do_not_return=False):

  assert mmseqs2_fn is not None
  if not do_not_filter:
    assert hhfilter_fn is not None

  if cov > 1: cov = cov/100

  mode = mode.split("_") if isinstance(mode,str) else mode
  seqs = [seqs] if isinstance(seqs,str) else seqs
  Ls = [len(seq) for seq in seqs]

  path = os.path.join(jobname,"msa")
  os.makedirs(path, exist_ok=True)
  if output_a3m is None:
    output_a3m = f"{jobname}/msa.a3m"

  msa_a3m = []
  if "paired" in mode and len(seqs) > 1:
    if verbose: print("getting paired MSA")
    out_paired = mmseqs2_fn(seqs, f"{path}/", use_pairing=True)
    for a3m_lines in out_paired:
      n = -1
      for line in a3m_lines.split("\n"):
        if len(line) > 0:
          if line.startswith(">"):
            n += 1
            if len(msa_a3m) < (n + 1):
              msa_a3m.append([])
          else:
            msa_a3m[n].append(line)

  if "unpaired" in mode or len(seqs) == 1:
    if verbose: print("getting unpaired MSA")
    out = mmseqs2_fn(seqs,f"{path}/")
    for n,a3m_lines in enumerate(out):
      for line in a3m_lines.split("\n"):
        if not line.startswith(">") and len(line) > 0:
          xs = ['-'*l for l in Ls]
          xs[n] = line.rstrip()
          msa_a3m.append(xs)

  if not do_not_return:
    N = len(msa_a3m)
    if verbose: print("parsing msas")
    msa,deletion_matrix = [],[]
    for n in range(N):
      msa_a3m[n] = "".join(msa_a3m[n])
      msa_vec = []
      deletion_vec = []
      deletion_count = 0
      for j in msa_a3m[n]:
        if j.islower():
          deletion_count += 1
        else:
          deletion_vec.append(deletion_count)
          deletion_count = 0
          msa_vec.append(j)
      msa.append(msa_vec)  
      deletion_matrix.append(deletion_vec)

  if do_not_filter:
    if not do_not_return:
      final_msa = msa
      final_deletion_matrix = deletion_matrix
    with open(output_a3m,"w") as handle:
      handle.write(f">X\n{''.join(seqs)}\n")    
      for n,seq in enumerate(msa_a3m):
        handle.write(f">{n}\n{''.join(seq)}\n")            
  
  else:
    if verbose: print("gathering info")
    se = [[sum(Ls[:i]), sum(Ls[:i+1]), Ls[i]] for i in range(len(Ls))]
    labels = {}
    for n in range(N):
      ok = True
      label = ""
      for (s,e,l) in se:
        num_gaps = msa[n][s:e].count("-")
        if num_gaps == l:
          label += "0"
        elif (1 - num_gaps/l) > cov:
          label += "1"
        else:
          ok = False
          break
      if ok:
        if label not in labels:
          labels[label] = []
        labels[label].append(n)
    
    if verbose: print("filtering sequences")
    Ns = {}
    for label,sub_Ns in labels.items():
      with open(f"{path}/{label}.a3m","w") as handle:
        
        # first sequence
        seq = ['-'*l for l in Ls]
        for m,x in enumerate(label):
          if x == "1": seq[m] = seqs[m]
        seq = "".join(seq)
        handle.write(f">X\n{seq}\n")
        
        for n in sub_Ns:
          handle.write(f">{n}\n{msa_a3m[n]}\n")
      
      hhfilter_fn(f"{path}/{label}.a3m", f"{path}/{label}.out.a3m", id=id, qid=qid)
      Ns[label] = []
      with open(f"{path}/{label}.out.a3m","r") as handle:
        for line in handle:
          if line.startswith(">") and not line.startswith(">X"):
            n = int(line[1:])
            Ns[label].append(n)

    if verbose: print("selecting final sequences")
    with open(output_a3m,"w") as handle:
      
      # add first sequence
      final_msa = [list("".join(seqs))]
      final_deletion_matrix = [[0]*len(final_msa[0])]
      handle.write(f">X\n{''.join(seqs)}\n")    

      counts = np.array([0] * len(Ls))
      label_to_counts = [[k,np.array(list(k)).astype(int)] for k in Ns.keys()]
      while len(final_msa) < max_msa:
        (best_k,best_v,best_sco) = (None,None,None)
        for k,v in label_to_counts:
          if len(Ns[k]) > 0:
            sco = (counts+v).std() - v.sum()
            if best_k is None or (sco < best_sco):
              (best_k,best_v,best_sco) = (k,v,sco)
        if best_k is None:
          break
        else:
          counts += best_v
          n = Ns[best_k].pop(0)
          final_msa.append(msa[n])
          final_deletion_matrix.append(deletion_matrix[n])
          handle.write(f">{n}\n{''.join(msa_a3m[n])}\n")            
  if not do_not_return:
    final_msa = np.array([[aa_order.get(res,20) for res in seq] for seq in final_msa])
    final_deletion_matrix = np.array(final_deletion_matrix)
    return final_msa, final_deletion_matrix