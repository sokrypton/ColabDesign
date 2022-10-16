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
          "length":np.array(len_set),
          "pos":np.asarray(pos_set)}

def rewire(length, order=None, loops=0, offset=0):
  '''
  Given a list of segment [length]s, move them around given an [offset], [order] and [loop] lengths.
  The [order] of the segments and the length of [loops] between segments can be controlled.
  '''
  seg_len = [length] if isinstance(length,int) else length
  num_seg = len(seg_len)

  # define order of segments
  if order is None: order = list(range(num_seg))
  assert len(order) == num_seg

  # define loop lengths between segments
  loop_len = ([loops] * (num_seg - 1)) if isinstance(loops, int) else loops
  assert len(loop_len) == num_seg - 1

  # get positions we want to restrain/constrain within hallucinated protein 
  l,new_pos = offset,[]
  for n,i in enumerate(np.argsort(order)):
    new_pos.append(l + np.arange(seg_len[i]))
    if n < num_seg - 1: l += seg_len[i] + loop_len[n] 

  return np.concatenate([new_pos[i] for i in order])