import jax
import jax.numpy as jnp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec 

from colabdesign.shared.protein import _np_kabsch
from colabdesign.shared.utils import update_dict, Key
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
    model.set_opt(num_models=1, num_recycles=0)
    model.set_opt(con=dict(num=1)) or set_opt({"con":{"num":1}}) or set_opt("con",num=1)
    model.set_opt(lr=1, set_defaults=True)
    '''
    ks = list(kwargs.keys())
    self.set_args(**{k:kwargs.pop(k) for k in ks if k in self._args})
        
    if kwargs.pop("set_defaults", False):
      update_dict(self._opt, *args, **kwargs)

    update_dict(self.opt, *args, **kwargs)

  def set_args(self, **kwargs):
    '''
    set [arg]uments
    '''
    for k in ["best_metric", "traj_iter", "shuffle_first"]:
      if k in kwargs: self._args[k] = kwargs.pop(k)
            
    if "recycle_mode" in kwargs:
      ok_recycle_mode_swap = ["average","sample","first","last"]
      if kwargs["recycle_mode"] in ok_recycle_mode_swap and self._args["recycle_mode"] in ok_recycle_mode_swap:
        self._args["recycle_mode"] = kwargs.pop("recycle_mode")
      else:
        print(f"ERROR: use {self.__class__.__name__}(recycle_mode=...) to set the recycle_mode")
    
    if "optimizer" in kwargs:
      self.set_optimizer(kwargs.pop("optimizer"),
        learning_rate=kwargs.pop("learning_rate", None))

    ks = list(kwargs.keys())
    if len(ks) > 0:
      print(f"ERROR: the following args were not set: {ks}")

  def get_loss(self, x="loss"):
    '''output the loss (for entire trajectory)'''
    return np.array([loss[x] for loss in self._tmp["log"]])

  def save_pdb(self, filename=None, get_best=True, renum_pdb=True, aux=None):
    '''
    save pdb coordinates (if filename provided, otherwise return as string)
    - set get_best=False, to get the last sampled sequence
    '''
    if aux is None:
      aux = self._tmp["best"]["aux"] if (get_best and "aux" in self._tmp["best"]) else self.aux
    aux = aux["all"]
    
    p = {k:aux[k] for k in ["aatype","residue_index","atom_positions","atom_mask"]}
    p["b_factors"] = 100 * p["atom_mask"] * aux["plddt"][...,None]

    def to_pdb_str(x, n=None):
      p_str = protein.to_pdb(protein.Protein(**x))
      p_str = "\n".join(p_str.splitlines()[1:-2])
      if renum_pdb: p_str = renum_pdb_str(p_str, self._lengths)
      if n is not None:
        p_str = f"MODEL{n:8}\n{p_str}\nENDMDL\n"
      return p_str

    p_str = ""
    for n in range(p["atom_positions"].shape[0]):
      p_str += to_pdb_str(jax.tree_map(lambda x:x[n],p), n+1)
    p_str += "END\n"
    
    if filename is None:
      return p_str
    else: 
      with open(filename, 'w') as f:
        f.write(p_str)

  #-------------------------------------
  # plotting functions
  #-------------------------------------
  def animate(self, s=0, e=None, dpi=100, get_best=True, traj=None, aux=None, color_by="plddt"):
    '''
    animate the trajectory
    - use [s]tart and [e]nd to define range to be animated
    - use dpi to specify the resolution of animation
    - color_by = ["plddt","chain","rainbow"]
    '''
    if aux is None:
      aux = self._tmp["best"]["aux"] if (get_best and "aux" in self._tmp["best"]) else self.aux
    aux = aux["all"]    
    if self.protocol in ["fixbb","binder"]:
      pos_ref = self._inputs["batch"]["all_atom_positions"][:,1].copy()
      pos_ref[(pos_ref == 0).any(-1)] = np.nan
    else:
      pos_ref = aux["atom_positions"][0,:,1,:]

    if traj is None: traj = self._tmp["traj"]
    sub_traj = {k:v[s:e] for k,v in traj.items()}

    align_xyz = self.protocol == "hallucination"
    return make_animation(**sub_traj, pos_ref=pos_ref, length=self._lengths,
                          color_by=color_by, align_xyz=align_xyz, dpi=dpi) 

  def plot_pdb(self, show_sidechains=False, show_mainchains=False,
               color="pLDDT", color_HP=False, size=(800,480), animate=False,
               get_best=True, aux=None, pdb_str=None):
    '''
    use py3Dmol to plot pdb coordinates
    - color=["pLDDT","chain","rainbow"]
    '''
    if pdb_str is None:
      pdb_str = self.save_pdb(get_best=get_best, aux=aux)
    view = show_pdb(pdb_str,
                    show_sidechains=show_sidechains,
                    show_mainchains=show_mainchains,
                    color=color,
                    Ls=self._lengths,
                    color_HP=color_HP,
                    size=size,
                    animate=animate)
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
      seqid = self.get_loss("seqid")
      ax1_.plot(seqid,color="green",label="seqid")
      # axes labels
      ax1.set_yscale("log")
      ticks = [0.25,0.5,1,2,4,8,16,32,64]
      ax1.set(xticks=[])
      ax1.set_yticks(ticks); ax1.set_yticklabels(ticks)
      ax1.set_ylabel("RMSD",color="black");ax1_.set_ylabel("seqid",color="green")
      ax1.set_ylim(0.25,64)
      ax1_.set_ylim(0,0.8)
      # extras
      ax2.plot(self.get_loss("soft"),color="yellow",label="soft")
      ax2.plot(self.get_loss("temp"),color="orange",label="temp")
      ax2.plot(self.get_loss("hard"),color="red",label="hard")
      ax2.set_ylim(-0.1,1.1)
      ax2.set_xlabel("iterations")
      ax2.legend(loc='center left')
    else:
      print("TODO")
    plt.show()

  def clear_best(self):
    self._tmp["best"] = {}

  def save_current_pdb(self, filename=None):
    '''save pdb coordinates (if filename provided, otherwise return as string)'''
    self.save_pdb(filename=filename, get_best=False)

  def plot_current_pdb(self, show_sidechains=False, show_mainchains=False,
    color="pLDDT", color_HP=False, size=(800,480), animate=False):
    '''use py3Dmol to plot pdb coordinates
    - color=["pLDDT","chain","rainbow"]
    '''
    self.plot_pdb(show_sidechains=show_sidechains, show_mainchains=show_mainchains, color=color,
      color_HP=color_HP, size=size, animate=animate, get_best=False)