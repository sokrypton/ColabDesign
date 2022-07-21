import jax
import jax.numpy as jnp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec 

from colabdesign.shared.protein import _np_kabsch
from colabdesign.shared.utils import update_dict
from colabdesign.shared.plot import plot_pseudo_3D, make_animation, show_pdb
from colabdesign.shared.protein import renum_pdb_str
from colabdesign.af.alphafold.common import protein

####################################################
# AF_UTILS - various utils (save, plot, etc)
####################################################
class _af_utils:  

  def set_opt(self, *args, **kwargs):
    '''
    set [opt]ions
    -------------------
    note: model.restart() resets the [opt]ions to their defaults
    use model.set_opt(..., set_defaults=True) 
    or model.restart(..., reset_opt=False) to avoid this
    -------------------    
    model.set_opt(models=1, recycles=0)
    model.set_opt(con=dict(num=1)) or set_opt({"con":{"num":1}})
    model.set_opt(lr=1, set_defaults=True)
    '''
    if "best_metric" in kwargs:
      self._args["best_metric"] = kwargs.pop("best_metric")

    if "optimizer" in kwargs:
      print("ERROR: use model.restart(optimizer=...) to set the optimizer")

    if "crop_mode" in kwargs:
      assert kwargs["crop_mode"] in ["slide","roll"]
      self._args["crop_mode"] = kwargs.pop("crop_mode")

    if "recycle_mode" in kwargs:
      if kwargs["recycle_mode"] in ["sample","last"] and self._args["recycle_mode"] in ["sample","last"]:
        self._args["recycle_mode"] = kwargs.pop("recycle_mode")
      else:
        print(f"ERROR: use {self.__class__.__name__}(recycle_mode=...) to set the recycle_mode")

    if kwargs.pop("set_defaults", False):
      update_dict(self._opt, *args, **kwargs)

    update_dict(self.opt, *args, **kwargs)
  
  def get_loss(self, x="loss"):
    '''output the loss (for entire trajectory)'''
    return np.array([float(loss[x]) for loss in self._traj["log"]])

  def save_pdb(self, filename=None, get_best=True, renum_pdb=True):
    '''
    save pdb coordinates (if filename provided, otherwise return as string)
    - set get_best=False, to get the last sampled sequence
    '''
    aux = self.aux if (self._best_aux is None or not get_best) else self._best_aux
    aux = jax.tree_map(np.asarray, aux)

    Ls = None    
    if self.protocol == "binder":
      Ls = [self._target_len, self._binder_len]      
    if self.protocol == "partial":
      Ls = [self._len]
    if self.protocol in ["hallucination","fixbb"]:
      Ls = [self._len] * self._copies
    
    p = {k:aux[k] for k in ["aatype","residue_index","atom_positions","atom_mask"]}        
    if p["aatype"].ndim == 2: p["aatype"] = p["aatype"].argmax(-1)

    b_factors = 100.0 * aux["plddt"][:,None] * p["atom_mask"]
    p = protein.Protein(**p,b_factors=b_factors)
    pdb_str = protein.to_pdb(p)
    if renum_pdb:
      pdb_str = renum_pdb_str(pdb_str, Ls)    
    if filename is None:
      return pdb_str, Ls
    else:
      with open(filename, 'w') as f: f.write(pdb_str)

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
    pos_ref = aux["atom_positions"][:,1,:]
    sub_traj = {k:v[s:e] for k,v in self._traj.items()}      
    if self.protocol == "hallucintion":
      length = [self._len] * self._copies
      return make_animation(**sub_traj, pos_ref=pos_ref, length=length, dpi=dpi)
    else:
      return make_animation(**sub_traj, pos_ref=pos_ref, align_xyz=False, dpi=dpi)  

  def plot_pdb(self, show_sidechains=False, show_mainchains=False,
               color="pLDDT", color_HP=False, size=(800,480), get_best=True):
    '''
    use py3Dmol to plot pdb coordinates
    - color=["pLDDT","chain","rainbow"]
    '''
    pdb_str, Ls = self.save_pdb(get_best=get_best)
    view = show_pdb(pdb_str,
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
    
    if self.protocol in ["fixbb","partial"] or (self.protocol == "binder" and self._args["redesign"]):
      if self.protocol == "partial" and self._args["use_sidechains"]:
        rmsd = self.get_loss("sc_rmsd")
      else:
        rmsd = self.get_loss("rmsd")
      for k in [0.5,1,2,4,8,16,32]:
        ax1.plot([0,len(rmsd)],[k,k],color="lightgrey")
      ax1.plot(rmsd,color="black")
      ax1_.plot(self.get_loss("seqid"),color="green",label="seqid")
      # axes labels
      ax1.set_yscale("log")
      ticks = [0.25,0.5,1,2,4,8,16,32,64]
      ax1.set(xticks=[])
      ax1.set_yticks(ticks); ax1.set_yticklabels(ticks)
      ax1.set_ylabel("RMSD",color="black");ax1_.set_ylabel("seqid",color="green")
      ax1.set_ylim(0.25,64)
      ax1_.set_ylim(0,0.4)
      # extras
      if "soft" in self._traj["log"][0]:
        ax2.plot(self.get_loss("soft"),color="yellow",label="soft")
      if "temp" in self._traj["log"][0]:
        ax2.plot(self.get_loss("temp"),color="orange",label="temp")
      if "hard" in self._traj["log"][0]:
        ax2.plot(self.get_loss("hard"),color="red",label="hard")
      ax2.set_ylim(-0.1,1.1)
      ax2.set_xlabel("iterations")
      ax2.legend(loc='center left')
    else:
      print("TODO")
    plt.show()