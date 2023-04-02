from string import ascii_uppercase, ascii_lowercase
alphabet_list = list(ascii_uppercase+ascii_lowercase)
import numpy as np

def sym_it(coords, center, cyclic_symmetry_axis, reflection_axis=None):

  def rotation_matrix(axis, theta):
      axis = axis / np.linalg.norm(axis)
      a = np.cos(theta / 2)
      b, c, d = -axis * np.sin(theta / 2)
      return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                      [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                      [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

  def align_axes(coords, source_axis, target_axis):
      rotation_axis = np.cross(source_axis, target_axis)
      rotation_angle = np.arccos(np.dot(source_axis, target_axis))
      rot_matrix = rotation_matrix(rotation_axis, rotation_angle)
      return np.dot(coords, rot_matrix)

  # Center the coordinates
  coords = coords - center

  # Align cyclic symmetry axis with Z-axis
  z_axis = np.array([0, 0, 1])
  coords = align_axes(coords, cyclic_symmetry_axis, z_axis)

  if reflection_axis is not None:
    # Align reflection axis with X-axis
    x_axis = np.array([1, 0, 0])
    coords = align_axes(coords, reflection_axis, x_axis)
  return coords


def fix_partial_contigs(contigs, parsed_pdb):
  INF = float("inf")

  # get unique chains
  chains = []
  for c, i in parsed_pdb["pdb_idx"]:
    if c not in chains: chains.append(c)

  # get observed positions and chains
  ok = []
  for contig in contigs:
    for x in contig.split("/"):
      if x[0].isalpha:
        C,x = x[0],x[1:]
        S,E = -INF,INF
        if x.startswith("-"):
          E = int(x[1:])
        elif x.endswith("-"):
          S = int(x[:-1])
        elif "-" in x:
          (S,E) = (int(y) for y in x.split("-"))
        elif x.isnumeric():
          S = E = int(x)      
        for c, i in parsed_pdb["pdb_idx"]:
          if c == C and i >= S and i <= E:
            if [c,i] not in ok: ok.append([c,i])

  # define new contigs
  new_contigs = []
  for C in chains:
    new_contig = []
    unseen = []
    seen = []
    for c,i in parsed_pdb["pdb_idx"]:
      if c == C:
        if [c,i] in ok:
          L = len(unseen)
          if L > 0:
            new_contig.append(f"{L}-{L}")
            unseen = []
          seen.append([c,i])
        else:
          L = len(seen)
          if L > 0:
            new_contig.append(f"{seen[0][0]}{seen[0][1]}-{seen[-1][1]}")
            seen = []
          unseen.append([c,i])
    L = len(unseen)
    if L > 0:
      new_contig.append(f"{L}-{L}")
    L = len(seen)
    if L > 0:
      new_contig.append(f"{seen[0][0]}{seen[0][1]}-{seen[-1][1]}")
    new_contigs.append("/".join(new_contig))

  return new_contigs

def fix_contigs(contigs,parsed_pdb):
  def fix_contig(contig):
    INF = float("inf")
    X = contig.split("/")
    Y = []
    for n,x in enumerate(X):
      if x[0].isalpha():
        C,x = x[0],x[1:]
        S,E = -INF,INF
        if x.startswith("-"):
          E = int(x[1:])
        elif x.endswith("-"):
          S = int(x[:-1])
        elif "-" in x:
          (S,E) = (int(y) for y in x.split("-"))
        elif x.isnumeric():
          S = E = int(x)      
        new_x = ""
        c_,i_ = None,0
        for c, i in parsed_pdb["pdb_idx"]:
          if c == C and i >= S and i <= E:
            if c_ is None:
              new_x = f"{c}{i}"
            else:
              if c != c_ or i != i_+1:
                new_x += f"-{i_}/{c}{i}"
            c_,i_ = c,i
        Y.append(new_x + f"-{i_}")
      elif "-" in x:
        # sample length
        s,e = x.split("-")
        m = np.random.randint(int(s),int(e)+1)
        Y.append(f"{m}-{m}")
      elif x.isnumeric() and x != "0":
        Y.append(f"{x}-{x}")
    return "/".join(Y)
  return [fix_contig(x) for x in contigs]

def fix_pdb(pdb_str, contigs):
  def get_range(contig):
    L_init = 1
    R = []
    sub_contigs = [x.split("-") for x in contig.split("/")]
    for n,(a,b) in enumerate(sub_contigs):
      if a[0].isalpha():
        if n > 0:
          pa,pb = sub_contigs[n-1]
          if pa[0].isalpha() and a[0] == pa[0]:
            L_init += int(a[1:]) - int(pb) - 1
        L = int(b)-int(a[1:]) + 1
      else:
        L = int(b)
      R += range(L_init,L_init+L)  
      L_init += L
    return R
  
  contig_ranges = [get_range(x) for x in contigs]
  R,C = [],[]
  for n,r in enumerate(contig_ranges):
    R += r
    C += [alphabet_list[n]] * len(r)
  
  pdb_out = []
  r_, c_,n = None, None, 0 
  for line in pdb_str.split("\n"):
    if line[:4] == "ATOM":
      c = line[21:22]
      r = int(line[22:22+5])
      if r_ is None: r_ = r
      if c_ is None: c_ = c
      if r != r_ or c != c_:
        n += 1
        r_,c_ = r,c
      pdb_out.append("%s%s%4i%s" % (line[:21],C[n],R[n],line[26:]))
    if line[:5] == "MODEL" or line[:3] == "TER" or line[:6] == "ENDMDL":
      pdb_out.append(line)
      r_, c_,n = None, None, 0 
  return "\n".join(pdb_out)