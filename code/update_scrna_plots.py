import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns
from plotting_functions import format_pvalue
import numpy as np
import pandas as pd

def leiden_cell_abundance_stripplots(adata_dict, leiden_col = 'leiden'):
    figs = []
    for tissue in adata_dict:    
        d = adata_dict[tissue].obs.value_counts(['sample', leiden_col]).unstack()
        d = (d.T/d.sum(axis=1)).T
        d['Genotype'] = d.index.str.split("_").str[1]
        fig = plt.figure(figsize = (12, 3))
    
        sns.stripplot(data = d.melt('Genotype').fillna(0), hue='Genotype', x=leiden_col, y='value',
                     palette = ['#FF6600', '#00CC99']
                     )
        pvals = scipy.stats.ttest_ind(d[d['Genotype']=='TIR1'].fillna(0).drop(['Genotype'], axis=1),
                          d[d['Genotype']=='WT'].fillna(0).drop(['Genotype'], axis=1),
                         )[1]
        pvals = scipy.stats.false_discovery_control(pvals)
        for c, i in enumerate(pvals):
            if i < .05:
                pval_str = format_pvalue(i)
                plt.text(c, 2e-1, pval_str, ha = 'center')
        plt.legend(bbox_to_anchor=(1, 1), loc='upper left', title='Genotype')
        plt.title(tissue)
        plt.yscale('log')
        plt.ylabel("Fraction of cells")
        figs.append(fig)
    return figs


def leiden_cell_abundance_stripplots_matched_t_test(adata_dict, leiden_col = 'leiden'):
    figs = []
    print("Using paired t-test")
    for tissue in adata_dict:    
        d = adata_dict[tissue].obs.value_counts(['sample', leiden_col]).unstack()
        d = (d.T/d.sum(axis=1)).T
        d['Genotype'] = d.index.str.split("_").str[1]
        fig = plt.figure(figsize = (12, 3))
    
        sns.stripplot(data = d.melt('Genotype').fillna(0), hue='Genotype', x=leiden_col, y='value',
                     palette = ['#FF6600', '#00CC99']
                     )
        
        values_1 = d[d['Genotype']=='TIR1'].fillna(0).drop(['Genotype'], axis=1)
        values_2 = d[d['Genotype']=='WT'].fillna(0).drop(['Genotype'], axis=1)
        # pvals = scipy.stats.ttest_ind(values_1, values_2)[1]
        pvals = scipy.stats.ttest_rel(values_1, values_2)[1]
        pvals = scipy.stats.false_discovery_control(pvals)
        for c, i in enumerate(pvals):
            if i < .05:
                pval_str = format_pvalue(i)
                plt.text(c, 2e-1, pval_str, ha = 'center')
        plt.legend(bbox_to_anchor=(1, 1), loc='upper left', title='Genotype')
        plt.title(tissue)
        plt.yscale('log')
        plt.ylabel("Fraction of cells")
        figs.append(fig)
    return figs


def n_deg_plot(deseq_padj_df, deseq_l2fc_df, pco=.05, column_order = None):
    n_degs = ((deseq_padj_df < pco) * np.sign(deseq_l2fc_df)).melt().value_counts().unstack().drop([0], axis=1)
    n_degs.columns = pd.Series(n_degs.columns).apply({-1 : 'Tir1-down', 1: 'Tir1-up'}.get)
    n_degs = n_degs.T
    n_degs = n_degs.loc[['Tir1-up', 'Tir1-down',]]
    if column_order is not None:
        n_degs = n_degs.loc[column_order]
    n_degs = n_degs.T
    n_degs.plot.bar(stacked=True, color=['tab:red', 'tab:blue'], linewidth=1, edgecolor='white')
    plt.xticks(rotation=0)
    plt.legend(bbox_to_anchor=(1, 1), title='# DEGs')
    plt.ylabel("# DEGs")
    plt.xlabel("Context")
    plt.title("# DEGs by context")
    return plt.gcf()    