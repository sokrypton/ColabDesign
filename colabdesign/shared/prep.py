import numpy as np
def prep_pos(pos, residue, chain):
  '''
  given input [pos]itions (a string of segment ranges seperated by comma,
  for example: "1,3-4,10-15"), return list of indices to constrain.
  '''
  residue_set = []
  chain_set = []
  len_set = []
  for idx in pos.split(","):
    i,j = idx.split("-") if "-" in idx else (idx, None)

    if i.isalpha() and j is None:
      residue_set += [None]
      chain_set += [i]
      len_set += [i]
    else:
      # if chain defined
      if i[0].isalpha():
        c,i = i[0], int(i[1:])
      else:
        c,i = chain[0],int(i)
      if j is None:
        j = i
      else:
        j = int(j[1:] if j[0].isalpha() else j)
      residue_set += list(range(i,j+1))
      chain_set += [c] * (j-i+1)
      len_set += [j-i+1]

  residue = np.asarray(residue)
  chain = np.asarray(chain)
  pos_set = []
  for i,c in zip(residue_set, chain_set):
    if i is None:
      idx = np.where(chain == c)[0]
      assert len(idx) > 0, f'ERROR: chain {c} not found'
      pos_set += [n for n in idx]
      len_set[len_set.index(c)] = len(idx)
    else:
      idx = np.where((residue == i) & (chain == c))[0]
      assert len(idx) == 1, f'ERROR: positions {i} and chain {c} not found'
      pos_set.append(idx[0])

  return {"residue":np.array(residue_set),
          "chain":np.array(chain_set),
          "seg_lengths":np.array(len_set),
          "pos":np.asarray(pos_set)}