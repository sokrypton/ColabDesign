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
from collections import OrderedDict, Counter

pymol_color_list = ["#33ff33","#00ffff","#ff33cc","#ffff00","#ff9999","#e5e5e5","#7f7fff","#ff7f00",
                    "#7fff7f","#199999","#ff007f","#ffdd5e","#8c3f99","#b2b2b2","#007fff","#c4b200",
                    "#8cb266","#00bfbf","#b27f7f","#fcd1a5","#ff7f7f","#ffbfdd","#7fffff","#ffff7f",
                    "#00ff7f","#337fcc","#d8337f","#bfff3f","#ff7fff","#d8d8ff","#3fffbf","#b78c4c",
                    "#339933","#66b2b2","#ba8c84","#84bf00","#b24c66","#7f7f7f","#3f3fa5","#a5512b"]

alphabet_list = list(ascii_uppercase+ascii_lowercase)
order_aa = {b: a for a, b in residue_constants.restype_order.items()}

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