# standard imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# logomaker import
import logomaker

def plot_with_logomaker(seq_vals,struct_vals,out,figsize=(20,3),title=None,mut_pos=None,editing_site_pos=None,xlim=None):
    f,axes=plt.subplots(nrows=2,ncols=1,dpi=80,figsize=figsize)
    seq_logo=logomaker.Logo(seq_vals,font_name='Arial Rounded MT Bold',color_scheme='classic',ax=axes[0])
    struct_logo=logomaker.Logo(struct_vals,font_name='Arial Rounded MT Bold',color_scheme={'S':'black','H':'black','I':'black','B':'black','E':'black','M':'black'},ax=axes[1])
    # style using Logo methods
    seq_logo.style_spines(visible=True)
    seq_logo.style_spines(spines=['left', 'bottom'], visible=True)
    seq_logo.style_xticks(rotation=90, fmt='%d', anchor=0, spacing=10)
    struct_logo.style_spines(visible=True)
    struct_logo.style_spines(spines=['left', 'bottom'], visible=True)
    struct_logo.style_xticks(rotation=90, fmt='%d', anchor=0, spacing=10)
    
    # style using Axes methods
    seq_logo.ax.set_ylabel("Sequence SHAP", labelpad=-1)
    seq_logo.ax.xaxis.set_tick_params(pad=-1)
    seq_logo.ax.plot(mut_pos, [-0.001]*len(mut_pos), '^r', markersize=15,alpha=0.4)
    seq_logo.ax.plot(editing_site_pos, -0.001, '^b', markersize=15,alpha=0.4)
    
    struct_logo.ax.set_ylabel("bpRNA Struct SHAP", labelpad=-1)
    struct_logo.ax.xaxis.set_tick_params(pad=-1)
    struct_logo.ax.plot(mut_pos, [-0.001]*len(mut_pos), '^r', markersize=15,alpha=0.4)
    struct_logo.ax.plot(editing_site_pos, -0.001, '^b', markersize=15,alpha=0.4)
    if xlim is not None:
        axes[0].set_xlim(xlim[0],xlim[1])
        axes[1].set_xlim(xlim[0],xlim[1])
    plt.suptitle(title)
    plt.subplots_adjust(hspace=0.4)
    plt.savefig(out,type='svg')
    

def plot_wrapper(seq,seq_shap,struct,struct_shap,out,figsize=(20,3),title=None,mut_pos=None, editing_site_pos=None,xlim=None):
    #make data frames of masked sequence & structure
    toplot_seq=pd.DataFrame(seq*seq_shap,columns=['A','C','G','U'])
    toplot_struct=pd.DataFrame(struct*struct_shap,columns=['S','H','B','I','E','M'])    
    #generate the plot
    plot_with_logomaker(toplot_seq,toplot_struct,out,figsize=figsize,title=title,mut_pos=mut_pos, editing_site_pos=editing_site_pos,xlim=xlim)
    
