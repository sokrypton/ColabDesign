import jax
import jax.numpy as jnp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec 

from colabdesign.shared.protein import _np_kabsch
from colabdesign.shared.utils import update_dict
from colabdesign.shared.plot import plot_pseudo_3D, make_animation, show_pdb
from colabdesign.af.alphafold.common import protein, residue_constants

ORDER_RESTYPE = {v: k for k, v in residue_constants.restype_order.items()}

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
    return ["".join([ORDER_RESTYPE[a] for a in s]) for s in x]
  
  def get_loss(self, x="loss"):
    '''output the loss (for entire trajectory)'''
    return np.array([float(loss[x]) for loss in self._traj["losses"]])

  def save_pdb(self, filename=None, get_best=True):
    '''
    save pdb coordinates (if filename provided, otherwise return as string)
    - set get_best=False, to get the last sampled sequence
    '''
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
      
    if self.protocol == "fixbb":
      pos_ref = self._batch["all_atom_positions"][:,1,:]
      return make_animation(**sub_traj, pos_ref=pos_ref, length=length, dpi=dpi)
    
    if self.protocol == "binder":
      pos_ref = aux["final_atom_positions"][:,1,:]
      return make_animation(**sub_traj, pos_ref=pos_ref, length=length, dpi=dpi)
    
    if self.protocol in ["hallucination","partial"]:
      pos_ref = aux["final_atom_positions"][:,1,:]
      return make_animation(**sub_traj, pos_ref=pos_ref, length=length, dpi=dpi)      
      

  def plot_pdb(self, show_sidechains=False, show_mainchains=False,
               color="pLDDT", color_HP=False, size=(800,480), get_best=True):
    '''
    use py3Dmol to plot pdb coordinates
    - color=["pLDDT","chain","rainbow"]
    '''
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