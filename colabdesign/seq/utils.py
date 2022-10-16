import os, string
import numpy as np
import jax
import jax.numpy as jnp

ALPHABET = list("ARNDCQEGHILKMFPSTWYV-")

def parse_fasta(filename, a3m=False, stop=100000):
  '''function to parse fasta file'''
  
  if a3m:
    # for a3m files the lowercase letters are removed
    # as these do not align to the query sequence
    rm_lc = str.maketrans(dict.fromkeys(string.ascii_lowercase))
    
  header, sequence = [],[]
  lines = open(filename, "r")
  for line in lines:
    line = line.rstrip()
    if len(line) > 0:
      if line[0] == ">":
        if len(header) == stop:
          break
        else:
          header.append(line[1:])
          sequence.append([])
      else:
        if a3m: line = line.translate(rm_lc)
        else: line = line.upper()
        sequence[-1].append(line)
  lines.close()
  sequence = [''.join(seq) for seq in sequence]
  
  return header, sequence
  
def mk_msa(seqs):
  '''one hot encode msa'''
  states = len(ALPHABET)  
  a2n = {a:n for n,a in enumerate(ALPHABET)}
  msa_ori = np.array([[a2n.get(aa, states-1) for aa in seq] for seq in seqs])
  return np.eye(states)[msa_ori]

def get_eff(msa, eff_cutoff=0.8):
  '''compute weight per sequence'''
  if msa.shape[0] > 10000:
    # loop one-to-all (to avoid memory issues)
    msa = msa.argmax(-1)
    def get_w(seq): return 1/((seq==msa).mean(-1) > eff_cutoff).sum()
    return jax.lax.scan(lambda _,x:(_,get_w(x)),None,msa,unroll=2)[1]
  else:
    # all-to-all
    msa_ident = jnp.tensordot(msa,msa,[[1,2],[1,2]])/msa.shape[1]
    return 1/(msa_ident >= eff_cutoff).sum(-1)   

def ar_mask(order, diag=True):
  '''compute autoregressive mask, given order of positions'''
  L = order.shape[0]
  r = order[::-1].argsort()
  tri = jnp.triu(jnp.ones((L,L)),k=not diag)
  return tri[r[None,:],r[:,None]]
