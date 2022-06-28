import jax
import jax.numpy as jnp
import numpy as np

from colabdesign.af.misc import jalview_color_list, _np_kabsch, order_restype
from colabdesign.af.alphafold.common import protein

# import matplotlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects
from matplotlib import animation
from matplotlib.gridspec import GridSpec 
from matplotlib import collections as mcoll

try:
  import py3Dmol
except:
  print("py3Dmol not installed")
  

from string import ascii_uppercase,ascii_lowercase
alphabet_list = list(ascii_uppercase+ascii_lowercase)
pymol_color_list = ["#33ff33","#00ffff","#ff33cc","#ffff00","#ff9999","#e5e5e5","#7f7fff","#ff7f00",
                    "#7fff7f","#199999","#ff007f","#ffdd5e","#8c3f99","#b2b2b2","#007fff","#c4b200",
                    "#8cb266","#00bfbf","#b27f7f","#fcd1a5","#ff7f7f","#ffbfdd","#7fffff","#ffff7f",
                    "#00ff7f","#337fcc","#d8337f","#bfff3f","#ff7fff","#d8d8ff","#3fffbf","#b78c4c",
                    "#339933","#66b2b2","#ba8c84","#84bf00","#b24c66","#7f7f7f","#3f3fa5","#a5512b"]
pymol_cmap = matplotlib.colors.ListedColormap(pymol_color_list)

# robust function for updating dictionary
def update_dict(D, *args, **kwargs):
  def set_dict(d, x):
    for k,v in x.items():
      if v is not None:
        if k in d:
          if isinstance(v,dict):
            set_dict(d[k], x[k])
          elif isinstance(d[k],dict):
            print(f"ERROR: '{k}' is a dictionary")
          elif isinstance(d[k],(int,float,bool,str)):
            d[k] = type(d[k])(v)
          else:
            d[k] = v
        else:
          print(f"ERROR: '{k}' not found in {list(d.keys())}")  
  for a in args:
    if isinstance(a, dict): set_dict(D, a)    
  set_dict(D, kwargs)

####################################################
# AF_UTILS - various utils (save, plot, etc)
####################################################
class _af_utils:
  def get_seqs(self, get_best=True):
    '''
    get sequences as strings
    - set get_best=False, to get the last sampled sequence
    '''
    aux = self.aux if (self._best_aux is None or not get_best) else self._best_aux
    aux = jax.tree_map(lambda x:np.asarray(x), aux)
    x = aux["seq"]["hard"].argmax(-1)
    return ["".join([order_restype[a] for a in s]) for s in x]
  
  def get_loss(self, x="loss"):
    '''output the loss (for entire trajectory)'''
    return np.array([float(loss[x]) for loss in self._traj["losses"]])

  def save_pdb(self, filename=None, get_best=True):
    '''
    save pdb coordinates (if filename provided, otherwise return as string)
    - set get_best=False, to get the last sampled sequence
    '''
    if self.use_struct:
      aux = self.aux if (self._best_aux is None or not get_best) else self._best_aux
      aux = jax.tree_map(lambda x:np.asarray(x), aux)
      aatype = aux["seq"]["hard"].argmax(-1)[0]
      if self.protocol == "binder":
        aatype_target = self._batch["aatype"][:self._target_len]
        aatype = np.concatenate([aatype_target,aatype])
      if self.protocol in ["fixbb","hallucination"] and self._copies > 1:
        aatype = np.concatenate([aatype] * self._copies)
      p = {"residue_index":self._inputs["residue_index"][0],
            "aatype":aatype,
            "atom_positions":aux["final_atom_positions"],
            "atom_mask":aux["final_atom_mask"]}
      b_factors = 100.0 * aux["plddt"][:,None] * p["atom_mask"]
      p = protein.Protein(**p,b_factors=b_factors)
      pdb_lines = protein.to_pdb(p)
      if filename is None:
        return pdb_lines
      else:
        with open(filename, 'w') as f: f.write(pdb_lines)
    else:
      print("ERROR: structure module disabled")
  #-------------------------------------
  # plotting functions
  #-------------------------------------
  def animate(self, s=0, e=None, dpi=100, get_best=True):
    '''
    animate the trajectory
    - use [s]tart and [e]nd to define range to be animated
    - use dpi to specify the resolution of animation
    '''
    aux = self.aux if (self._best_aux is None or not get_best) else self._best_aux
    
    sub_traj = {k:v[s:e] for k,v in self._traj.items()}
    length = None    
    if self.protocol == "binder":
      length = [self._target_len, self._binder_len]      
    if self.protocol == "partial":
      length = [self._len]
    if self.protocol in ["hallucination","fixbb"]:
      length = [self._len] * self._copies
      
    if self.use_struct:      
      
      if self.protocol == "fixbb":
        pos_ref = self._batch["all_atom_positions"][:,1,:]
        return make_animation(**sub_traj, pos_ref=pos_ref, length=length, dpi=dpi)
      
      if self.protocol == "binder":
        pos_ref = aux["final_atom_positions"][:,1,:]
        return make_animation(**sub_traj, pos_ref=pos_ref, length=length, dpi=dpi)
      
      if self.protocol in ["hallucination","partial"]:
        pos_ref = aux["final_atom_positions"][:,1,:]
        return make_animation(**sub_traj, pos_ref=pos_ref, length=length, dpi=dpi)      
      
    else:
      return make_animation(**sub_traj, length=length, dpi=dpi)

  def plot_pdb(self, show_sidechains=False, show_mainchains=False,
               color="pLDDT", color_HP=False, size=(800,480), get_best=True):
    '''
    use py3Dmol to plot pdb coordinates
    - color=["pLDDT","chain","rainbow"]
    '''
    if self.use_struct:
      Ls = None    
      if self.protocol == "binder":
        Ls = [self._target_len, self._binder_len]      
      if self.protocol == "partial":
        Ls = [self._len]
      if self.protocol in ["hallucination","fixbb"]:
        Ls = [self._len] * self._copies
      view = show_pdb(self.save_pdb(get_best=get_best),
                      show_sidechains=show_sidechains,
                      show_mainchains=show_mainchains,
                      color=color,
                      Ls=Ls,
                      color_HP=color_HP,
                      size=size)
      view.show()
    else:
      print("ERROR: structure module disabled")
  
  def plot_traj(self, dpi=100):
    fig = plt.figure(figsize=(5,5), dpi=dpi)
    gs = GridSpec(4,1, figure=fig)
    ax1 = fig.add_subplot(gs[:3,:])
    ax2 = fig.add_subplot(gs[3:,:])
    ax1_ = ax1.twinx()
    
    if self.protocol == "fixbb" or (self.protocol == "binder" and self._redesign):
      rmsd = self.get_loss("rmsd")
      for k in [0.5,1,2,4,8,16,32]:
        ax1.plot([0,len(rmsd)],[k,k],color="lightgrey")
      ax1.plot(rmsd,color="black")
      ax1_.plot(self.get_loss("seqid"),color="green",label="seqid")
      # axes labels
      ax1.set_yscale("log")
      ticks = [0.25,0.5,1,2,4,8,16,32,64]
      ax1.set(xticks=[])
      ax1.set_yticks(ticks);ax1.set_yticklabels(ticks)
      ax1.set_ylabel("RMSD",color="black");ax1_.set_ylabel("seqid",color="green")
      ax1.set_ylim(0.25,64)
      ax1_.set_ylim(0,0.4)
      # extras
      if "soft" in self._traj["losses"][0]:
        ax2.plot(self.get_loss("soft"),color="yellow",label="soft")
      if "temp" in self._traj["losses"][0]:
        ax2.plot(self.get_loss("temp"),color="orange",label="temp")
      if "hard" in self._traj["losses"][0]:
        ax2.plot(self.get_loss("hard"),color="red",label="hard")
      ax2.set_ylim(-0.1,1.1)
      ax2.set_xlabel("iterations")
      ax2.legend(loc='center left')
    else:
      print("TODO")
    plt.show()

def plot_pseudo_3D(xyz, c=None, ax=None, chainbreak=5,
                   cmap="gist_rainbow", line_w=2.0,
                   cmin=None, cmax=None, zmin=None, zmax=None):

  def rescale(a, amin=None, amax=None):
    a = np.copy(a)
    if amin is None: amin = a.min()
    if amax is None: amax = a.max()
    a[a < amin] = amin
    a[a > amax] = amax
    return (a - amin)/(amax - amin)

  # make segments
  xyz = np.asarray(xyz)
  seg = np.concatenate([xyz[:-1,None,:],xyz[1:,None,:]],axis=-2)
  seg_xy = seg[...,:2]
  seg_z = seg[...,2].mean(-1)
  ord = seg_z.argsort()

  # set colors
  if c is None: c = np.arange(len(seg))[::-1]
  else: c = (c[1:] + c[:-1])/2
  c = rescale(c,cmin,cmax)  

  if isinstance(cmap, str):
    if cmap == "gist_rainbow": c *= 0.75
    colors = matplotlib.cm.get_cmap(cmap)(c)
  else:
    colors = cmap(c)
  
  if chainbreak is not None:
    dist = np.linalg.norm(xyz[:-1] - xyz[1:], axis=-1)
    colors[...,3] = (dist < chainbreak).astype(np.float)

  # add shade/tint based on z-dimension
  z = rescale(seg_z,zmin,zmax)[:,None]
  tint, shade = z/3, (z+2)/3
  colors[:,:3] = colors[:,:3] + (1 - colors[:,:3]) * tint
  colors[:,:3] = colors[:,:3] * shade

  set_lim = False
  if ax is None:
    fig, ax = plt.subplots()
    fig.set_figwidth(5)
    fig.set_figheight(5)
    set_lim = True
  else:
    fig = ax.get_figure()
    if ax.get_xlim() == (0,1):
      set_lim = True
      
  if set_lim:
    xy_min = xyz[:,:2].min() - line_w
    xy_max = xyz[:,:2].max() + line_w
    ax.set_xlim(xy_min,xy_max)
    ax.set_ylim(xy_min,xy_max)

  ax.set_aspect('equal')
    
  # determine linewidths
  width = fig.bbox_inches.width * ax.get_position().width
  linewidths = line_w * 72 * width / np.diff(ax.get_xlim())

  lines = mcoll.LineCollection(seg_xy[ord], colors=colors[ord], linewidths=linewidths,
                               path_effects=[matplotlib.patheffects.Stroke(capstyle="round")])
  
  return ax.add_collection(lines)

def plot_ticks(ax, Ls, Ln=None, add_yticks=False):
  if Ln is None: Ln = sum(Ls)
  L_prev = 0
  for L_i in Ls[:-1]:
    L = L_prev + L_i
    L_prev += L_i
    ax.plot([0,Ln],[L,L],color="black")
    ax.plot([L,L],[0,Ln],color="black")
  
  if add_yticks:
    ticks = np.cumsum([0]+Ls)
    ticks = (ticks[1:] + ticks[:-1])/2
    ax.yticks(ticks,alphabet_list[:len(ticks)])
  
def make_animation(seq, con=None, xyz=None, plddt=None, pae=None,
                   losses=None, pos_ref=None, line_w=2.0,
                   dpi=100, interval=60, color_msa="Taylor",
                   length=None, **kwargs):

  def align(P, Q, P_trim=None):
    if P_trim is None: P_trim = P
    p_trim = P_trim - P_trim.mean(0,keepdims=True)
    p = P - P_trim.mean(0,keepdims=True)
    q = Q - Q.mean(0,keepdims=True)
    return p @ _np_kabsch(p_trim, q, use_jax=False)

  if xyz is not None:
    # compute reference position
    if pos_ref is None: pos_ref = xyz[-1]
    

    # align to reference
    if length is None:
      L = len(pos_ref)
    elif isinstance(length, list):
      L = length[0]
    else:
      L = length
      
    pos_ref_trim = pos_ref[:L]
    # align to reference position
    new_positions = []
    for i in range(len(xyz)):
      new_positions.append(align(xyz[i],pos_ref_trim,xyz[i][:L]))
    pos = np.asarray(new_positions)

    # rotate for best view
    pos_mean = np.concatenate(pos,0)
    m = pos_mean.mean(0)
    rot_mtx = _np_kabsch(pos_mean - m, pos_mean - m, return_v=True, use_jax=False)
    pos = (pos - m) @ rot_mtx + m
    pos_ref_full = ((pos_ref - pos_ref_trim.mean(0)) - m) @ rot_mtx + m

  # initialize figure
  if pae is not None and len(pae) == 0: pae = None
  fig = plt.figure()
  gs = GridSpec(4,3, figure=fig)
  if pae is not None:
    ax1, ax2, ax3 = fig.add_subplot(gs[:3,:2]), fig.add_subplot(gs[3:,:]), fig.add_subplot(gs[:3,2:])
  else:
    ax1, ax2 = fig.add_subplot(gs[:3,:]), fig.add_subplot(gs[3:,:])

  fig.subplots_adjust(top=0.95,bottom=0.1,right=0.95,left=0.05,hspace=0,wspace=0)
  fig.set_figwidth(8); fig.set_figheight(6); fig.set_dpi(dpi)
  ax2.set_xlabel("positions"); ax2.set_yticks([])
  if seq[0].shape[0] > 1: ax2.set_ylabel("sequences")
  else: ax2.set_ylabel("amino acids")

  if xyz is None:
    ax1.set_title("predicted contact map")
  else:
    ax1.set_title("N→C") if plddt is None else ax1.set_title("pLDDT")
  if pae is not None:
    ax3.set_title("pAE")
    ax3.set_xticks([])
    ax3.set_yticks([])

  # set bounderies
  if xyz is not None:
    x_min,y_min,z_min = np.minimum(np.mean(pos.min(1),0),pos_ref_full.min(0)) - 5
    x_max,y_max,z_max = np.maximum(np.mean(pos.max(1),0),pos_ref_full.max(0)) + 5

    x_pad = ((y_max - y_min) * 2 - (x_max - x_min)) / 2
    y_pad = ((x_max - x_min) / 2 - (y_max - y_min)) / 2
    if x_pad > 0:
      x_min -= x_pad
      x_max += x_pad
    else:
      y_min -= y_pad
      y_max += y_pad

    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
  ax1.set_xticks([])
  ax1.set_yticks([])

  # get animation frames
  ims = []
  for k in range(len(seq)):
    ims.append([])
    if xyz is not None:
      if plddt is None:
        ims[-1].append(plot_pseudo_3D(pos[k], ax=ax1, line_w=line_w, zmin=z_min, zmax=z_max))
      else:
        ims[-1].append(plot_pseudo_3D(pos[k], c=plddt[k], cmin=0.5, cmax=0.9, ax=ax1, line_w=line_w, zmin=z_min, zmax=z_max))
    else:
      L = con[k].shape[0]
      ims[-1].append(ax1.imshow(con[k], animated=True, cmap="Greys",vmin=0, vmax=1, extent=(0, L, L, 0)))

    if seq[k].shape[0] == 1:
      ims[-1].append(ax2.imshow(seq[k][0].T, animated=True, cmap="bwr_r",vmin=-1, vmax=1))
    else:
      cmap = matplotlib.colors.ListedColormap(jalview_color_list[color_msa])
      vmax = len(jalview_color_list[color_msa]) - 1
      ims[-1].append(ax2.imshow(seq[k].argmax(-1), animated=True, cmap=cmap, vmin=0, vmax=vmax, interpolation="none"))
    
    if pae is not None:
      L = pae[k].shape[0]
      ims[-1].append(ax3.imshow(pae[k], animated=True, cmap="bwr",vmin=0, vmax=30, extent=(0, L, L, 0)))

  # add lines
  if length is not None:
    Ls = length if isinstance(length, list) else [length,None]
    if con is not None:
      plot_ticks(ax1, Ls, con[0].shape[0])
    if pae is not None:
      plot_ticks(ax3, Ls, pae[0].shape[0])

  # make animation!
  ani = animation.ArtistAnimation(fig, ims, blit=True, interval=interval)
  plt.close()
  return ani.to_html5_video()

def clear_mem():
  backend = jax.lib.xla_bridge.get_backend()
  for buf in backend.live_buffers(): buf.delete()

def renum_pdb_str(pdb_str, Ls=None):
  if Ls is not None:
    L_init = 0
    new_chain = {}
    for L,c in zip(Ls, alphabet_list):
      new_chain.update({i:c for i in range(L_init,L_init+L)})
      L_init += L  

  n,pdb_out = 1,[]
  resnum_,chain_ = None,None
  for line in pdb_str.split("\n"):
    if line[:4] == "ATOM":
      chain = line[21:22]
      resnum = int(line[22:22+5])
      if resnum_ is None: resnum_ = resnum
      if chain_ is None: chain_ = chain
      if resnum != resnum_ or chain != chain_:
        resnum_,chain_ = resnum,chain
        n += 1
      if Ls is None: pdb_out.append("%s%4i%s" % (line[:22],n,line[26:]))
      else: pdb_out.append("%s%s%4i%s" % (line[:21],new_chain[n-1],n,line[26:]))        
  return "\n".join(pdb_out)
    
def show_pdb(pdb_str, show_sidechains=False, show_mainchains=False,
             color="pLDDT", chains=None, Ls=None, vmin=50, vmax=90,
             color_HP=False, size=(800,480), hbondCutoff=4.0):
  
  if chains is None:
    chains = 1 if Ls is None else len(Ls)

  view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js', width=size[0], height=size[1])
  pdb_str_renum = renum_pdb_str(pdb_str, Ls)
  view.addModel(pdb_str_renum,'pdb',{'hbondCutoff':hbondCutoff})
  if color == "pLDDT":
    view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':vmin,'max':vmax}}})
  elif color == "rainbow":
    view.setStyle({'cartoon': {'color':'spectrum'}})
  elif color == "chain":
    for n,chain,color in zip(range(chains),alphabet_list,pymol_color_list):
       view.setStyle({'chain':chain},{'cartoon': {'color':color}})
  if show_sidechains:
    BB = ['C','O','N']
    HP = ["ALA","GLY","VAL","ILE","LEU","PHE","MET","PRO","TRP","CYS","TYR"]
    if color_HP:
      view.addStyle({'and':[{'resn':HP},{'atom':BB,'invert':True}]},
                    {'stick':{'colorscheme':"yellowCarbon",'radius':0.3}})
      view.addStyle({'and':[{'resn':HP,'invert':True},{'atom':BB,'invert':True}]},
                    {'stick':{'colorscheme':"whiteCarbon",'radius':0.3}})
      view.addStyle({'and':[{'resn':"GLY"},{'atom':'CA'}]},
                    {'sphere':{'colorscheme':"yellowCarbon",'radius':0.3}})
      view.addStyle({'and':[{'resn':"PRO"},{'atom':['C','O'],'invert':True}]},
                    {'stick':{'colorscheme':"yellowCarbon",'radius':0.3}})
    else:
      view.addStyle({'and':[{'resn':["GLY","PRO"],'invert':True},{'atom':BB,'invert':True}]},
                    {'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
      view.addStyle({'and':[{'resn':"GLY"},{'atom':'CA'}]},
                    {'sphere':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
      view.addStyle({'and':[{'resn':"PRO"},{'atom':['C','O'],'invert':True}]},
                    {'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})  
  if show_mainchains:
    BB = ['C','O','N','CA']
    view.addStyle({'atom':BB},{'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
  view.zoomTo()
  return view
