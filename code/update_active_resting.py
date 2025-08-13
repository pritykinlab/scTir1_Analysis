from plotting_functions import init_subplots_exact
import scanpy as sc
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def make_active_resting_cutoffs(dividing_cutoff, active_cutoff, resting_cutoff, not_active_cutoff, adata):    
    adata.obs['cell_is_dividing'] = (adata.obs['G2M_score'] > dividing_cutoff).astype(int)

    adata.obs['cell_is_active'] = ((adata.obs['up_aTreg'] > active_cutoff)
                                   & ~(adata.obs['cell_is_dividing'])
                                  ).astype(int)
    
    fig, axs = init_subplots_exact(3, 1, fgsz=(4, 2), dpi = 50)
    plt.sca(axs[0])
    sns.histplot(adata.obs['G2M_score'], color='red', element='step', fill=False, label='G2M', stat='density')
    sns.histplot(adata.obs['S_score'], color='blue', element='step', fill=False, label='S', stat='density')
    plt.axvline(dividing_cutoff)
    
    plt.ylim([0, 3])
    plt.xlim(0, 1)

    plt.axvline(active_cutoff)
    plt.axvline(not_active_cutoff, color='red')

    
    plt.sca(axs[1])
    sns.histplot(adata.obs['up_aTreg'])
    plt.axvline(active_cutoff)
    plt.axvline(not_active_cutoff, color='red')

    adata.obs['cell_is_resting'] = ((adata.obs['up_rTreg'] > resting_cutoff) 
                                  & ((adata.obs['up_aTreg'] < not_active_cutoff) )
                                  & (~adata.obs['cell_is_active'])
                                  & ~(adata.obs['cell_is_dividing'])
                                   ).astype(int)
    plt.sca(axs[2])
    sns.histplot(adata.obs['up_rTreg'])
    plt.axvline(resting_cutoff)
    fig = sc.pl.umap(adata, color=[ 'up_aTreg', 'up_rTreg', 'G2M_score', 'cell_is_resting', 'cell_is_active', 'cell_is_dividing'],
               title=['aTreg Score', 'rTreg Score', 'G2M Score', 'Resting Treg', 'Active Treg', 'Dividing Treg'], 
               cmap='coolwarm', ncols = 3, vmax = 1, vmin = -1, return_fig=True, )
       
    return fig