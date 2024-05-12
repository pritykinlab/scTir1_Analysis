# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from plot_bulk_RNA import *
import sys
from foxp3_pileup_plots import plot_foxp3_peaks_near_degs_as_bar


arr = np.asarray
def make_gene_modules_on_all_genes(lfc_df_all, pval_df_all, basemean_df_all, 
                                   all_lfc_cols, all_pval_cols, 
                                   col_to_name, annot_cols,
                                   n_clust = 16
                                   ):
    basemean_co = 200
    genes_to_plot = arr(['Neb',  'Hspa1a', 'Lamc1', 'Snx33', 'Prss12', 'Il21', 'Id2', 'Itga5', 'Cd7', 'Dusp2', 'Socs1', 
                        'Izumo1r', 'Tnfrsf18', 'Gata3', 'Irf1', 'Tnfrsf4', 'Ly6e', 'Jun', 'Ifit2', 'Cd6', 'Lcp1', 
                        'Il6st', 'Tgfbr2', 'Fos', 'Dusp1', 'Tcf7', 'Tmem9b', 'Sox4', 'Tnfsf8', 'Dusp6', 'Pde3b', 
                        'Sema4a', 'Il12rb1', 'Coro2a', 'Setdb1', 'Il2ra', 'Il2rb', 'Foxj2', 'Tnfrsf1b', 'Lrrc32',])
    all_cluster_df_list = []
    all_cluster_df_color_dict = {}
    for key, tmp_annot_cols in all_lfc_cols.items():
        pval_cols = all_pval_cols[key]
        degs = pval_df_all.index[(pval_df_all[pval_cols] < 1).any(axis=1) 
                            & (basemean_df_all[pval_cols] > basemean_co).any(axis=1)]
        
        lfc_cols = all_lfc_cols[key]
        mat = lfc_df_all.loc[degs, lfc_cols]
        mat = mat.loc[~mat.isna().all(axis=1)]
        mat = mat.fillna(0).copy()
        mat.columns = [col_to_name.get(x, x) for x in mat.columns]
        mat = mat
        # mat = mat.clip(-2, 2)
        kmeans = KMeans(n_clusters=n_clust)
        kmeans.fit(mat)
        cluster2 = kmeans.labels_    
        o2 = np.argsort(cluster2)
        
        z = pd.concat([mat, pd.Series(cluster2, index=mat.index, name='Cluster')], axis=1).groupby("Cluster").mean().mean(axis=1).sort_values().index
        cluster2new = cluster2.copy()
        for c, i in enumerate(z):
            cluster2new[cluster2==i] = c
        cluster2 = cluster2new.copy()
        z = pd.concat([mat, pd.Series(cluster2, index=mat.index, name='Cluster')], axis=1).groupby("Cluster").mean().mean(axis=1).sort_values().index
        zvals = pd.concat([mat, pd.Series(cluster2, index=mat.index, name='Cluster')], axis=1).groupby("Cluster").mean().mean(axis=1).sort_values().values
        
        cluster_color_dict = dict(zip(z, sns.color_palette('coolwarm', as_cmap=True)((zvals + 2)/4)))
        cluster_quant_dict = dict(zip(z, zvals))

        all_cluster_df_color_dict[key] = cluster_color_dict
        
        fig, cluster_df = plot_bulk_rna_heatmap_and_clusters(mat, mat, lfc_df_all, pval_df_all, cluster2, 
                                                annot_cols, o2, genes_to_plot[np.isin(genes_to_plot, mat.index)],
                                                            cluster_color_dict,  cluster_quant_dict,
                                                            cluster_max = 2)
        plt.grid(False)
        plt.tight_layout()
        fig.axes[0].grid(False)    
        fig.axes[0].set_title(key)
        cluster2 = pd.Series(cluster2, name = key, index=mat.index)
        all_cluster_df_list.append(cluster2)
        fig.savefig(f'./FINAL_FIGURES/gene_modules/{key.replace(" ", "_")}.pdf', bbox_inches='tight')
    all_cluster_df = pd.concat(all_cluster_df_list, axis=1)
    return all_cluster_df, all_cluster_df_color_dict

from aux_functions import rename_clusts_by_clusts_in_order
def make_gene_modules_on_all_genes_wald(wald_df_all, pval_df_all, basemean_df_all, 
                                   all_lfc_cols, 
                                   col_to_name, annot_cols,
                                   n_clust = 16, basemean_co = 200, show=True):
    genes_to_plot = arr(['Neb',  'Hspa1a', 'Lamc1', 'Snx33', 'Prss12', 'Il21', 'Id2', 'Itga5', 'Cd7', 'Dusp2', 'Socs1', 
                        'Izumo1r', 'Tnfrsf18', 'Gata3', 'Irf1', 'Tnfrsf4', 'Ly6e', 'Jun', 'Ifit2', 'Cd6', 'Lcp1', 
                        'Il6st', 'Tgfbr2', 'Fos', 'Dusp1', 'Tcf7', 'Tmem9b', 'Sox4', 'Tnfsf8', 'Dusp6', 'Pde3b', 
                        'Sema4a', 'Il12rb1', 'Coro2a', 'Setdb1', 'Il2ra', 'Il2rb', 'Foxj2', 'Tnfrsf1b', 'Lrrc32',])
    all_cluster_df_list = []
    all_cluster_df_color_dict = {}
    all_cluster_quant_dict = {}
    for key, lfc_cols in all_lfc_cols.items():
        degs = basemean_df_all.index[(basemean_df_all[lfc_cols] > basemean_co).any(axis=1)]
        mat = wald_df_all.loc[degs, lfc_cols]
        mat = mat.loc[~mat.isna().all(axis=1)]
        mat = mat.dropna().copy()
        mat.columns = [col_to_name.get(x, x) for x in mat.columns]

        kmeans = KMeans(n_clusters=n_clust, random_state=0, n_init = 10, max_iter=1000)
        kmeans.fit(mat)
        cluster2 = kmeans.labels_    
        submat = mat.copy()
        submat['Cluster'] = cluster2
        # o2 = np.argsort(-cluster2)

        o2 = submat.groupby('Cluster').mean().mean(axis=1).sort_values().index
        cluster3 = rename_clusts_by_clusts_in_order(cluster2, submat.groupby('Cluster').mean().mean(axis=1).sort_values().index)
        o2 = np.argsort(cluster3)
        cluster2 = cluster3
        z = pd.concat([mat, pd.Series(cluster2, index=mat.index, name='Cluster')], axis=1).groupby("Cluster").mean().mean(axis=1).sort_values().index
        cluster2new = cluster2.copy()
        for c, i in enumerate(z):
            cluster2new[cluster2==i] = c
        cluster2 = cluster2new.copy()
        z = pd.concat([mat, pd.Series(cluster2, index=mat.index, name='Cluster')], axis=1).groupby("Cluster").mean().mean(axis=1).sort_values().index
        zvals = pd.concat([mat, pd.Series(cluster2, index=mat.index, name='Cluster')], axis=1).groupby("Cluster").mean().mean(axis=1).sort_values().values
        
        cluster_color_dict = dict(zip(z, sns.color_palette('coolwarm', as_cmap=True)((zvals + 4)/8)))
        cluster_quant_dict = dict(zip(z, zvals))

        all_cluster_quant_dict[key] = cluster_quant_dict
        all_cluster_df_color_dict[key] = cluster_color_dict

        fig, cluster_df = plot_bulk_rna_heatmap_and_clusters(mat, mat, wald_df_all, pval_df_all, cluster2, 
                                                annot_cols, o2, genes_to_plot[np.isin(genes_to_plot, mat.index)],
                                                            cluster_color_dict,  cluster_quant_dict,
                                                            cluster_max = 4,
                                                            cbar_label=r'Wald')
        plt.grid(False)
        fig.axes[0].grid(False)    
        fig.axes[0].set_title(key)
        cluster2 = pd.Series(cluster2, name = key, index=mat.index)
        all_cluster_df_list.append(cluster2)
        fig.savefig(f'./FINAL_FIGURES/gene_modules/{key.replace(" ", "_")}.pdf', bbox_inches='tight')
        if show==False:
            plt.close(fig)
    all_cluster_df = pd.concat(all_cluster_df_list, axis=1)
    return all_cluster_df, all_cluster_df_color_dict, all_cluster_quant_dict

def plot_gene_module_barplots(all_cluster_df, col_to_name, foxp3_count_set, foxp3_counts, all_cluster_df_color_dict, ymax=3,
                              prefix='gene_modules', title = '# Foxp3 peaks per gene'):
    n = len(all_cluster_df.columns)
    foxp3_count_set = foxp3_count_set & set(all_cluster_df.index)
    for c, col in enumerate(all_cluster_df.columns):
        fig, axs = init_subplots(1, 1, fgsz=(8, 3), space=.3)
        plt.sca(axs[0])
        us = all_cluster_df[col].dropna().astype(int).unique()
        ns = []
        for u in sorted(us):
            idx = all_cluster_df[col] == u
            degs = idx.index[idx]
            degs = degs[degs.isin(foxp3_counts.index)]
            if len(degs) < 3:
                continue

            v = plot_foxp3_peaks_near_degs_as_bar(degs, col_to_name, u, foxp3_count_set, foxp3_counts,
                                                xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding',
                                                color = all_cluster_df_color_dict[col][int(u)],
                                            linewidth = 3
                                )
            ns.append(len(degs))
        v = plot_foxp3_peaks_near_degs_as_bar(degs, col_to_name, u, foxp3_count_set, foxp3_counts,
                                            xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding',
                                            color = all_cluster_df_color_dict[col][int(u)],
                                        linewidth = 3, plot_baseline=True
                            )
        plt.xticks([-1] + list(range(0, len(us))))
        plt.ylim([0, ymax])
        plt.gca().set_xticklabels(['Baseline'] + list(range(0, len(us))))
        plt.title(f'{title} \n' + col)
        plt.ylabel("peaks per gene")
        plt.legend(bbox_to_anchor=(1.1, 1), ncol=2, frameon=False,
                   loc='upper left')
        
        xs = np.arange(len(ns))
        plt.twinx()
        plt.scatter(xs, ns, color='black', label='# Genes') 
        plt.grid(False)
        plt.ylim([0, 2500])
    #     plt.yticks([])
        plt.legend(loc='upper right', frameon=True)
        fig.savefig(f'./FINAL_FIGURES/gene_modules_foxp3_enrichment/{prefix}_{col.replace(" ", "_")}.pdf', bbox_inches='tight')    
    

from plotting_functions import add_xaxis_labels
from plotting_functions import init_subplots_exact
from plotting_functions import *


def plot_genemodules_on_reintro_data(cols, all_cluster_df, basemean_df_all, lfc_df_all, col_to_name, all_cluster_df_color_dict,
        basemean_co = 200):
    n = len(cols)
    fig, axs = init_subplots_exact(n, 1, fgsz=(100*mm, 100*mm), dpi = 100, as_list=True)
    for c, (wald_key, reintro_cols, reintro_cond) in enumerate(cols):
        plt.sca(axs[c])
        z = all_cluster_df[wald_key].dropna()
        for u in np.unique(z):
            genes = z.index[z == u]
            genes = genes[genes.isin(basemean_df_all.index)]
            genes = genes[(basemean_df_all.loc[genes, reintro_cols] > basemean_co).all(axis=1)]
            if len(genes) < 15:
                continue
            xs = [x.split('.')[1] for x in reintro_cols]
            data = lfc_df_all.loc[genes, reintro_cols].melt()#.mean(axis=0)
            data['variable'] = data['variable'].apply(lambda x: col_to_name.get(x).split(" ")[1])
            sns.lineplot(data=data, x= 'variable', y='value',
                    c = all_cluster_df_color_dict[wald_key][int(u)], label=u,legend=False,
                    )
        plt.title(wald_key.split("_vs_")[0])
        plt.ylim([-2, 2])
        plt.yticks([-2, -1, 0, 1, 2])
        plt.ylabel("LFC")
        plt.xlabel(f"Day, {reintro_cond} Reintro")
    fig.savefig('./FINAL_FIGURES/modules_on_bulk_reintro/lineplot.pdf', bbox_inches='tight')
        
        

import pandas as pd
def plot_peak_counts_by_wald_stat(cols, wald_df_all, basemean_df_all, foxp3_counts, col_to_name, ymax=3, factor='Foxp3', 
                                  n = 20):
    basemean_co = 50
    
    colors = ['red', 'blue', 'steelblue', 'orange', 'green', 'lightgreen']
    for c, col in enumerate(cols):
        idx = wald_df_all.index.isin(foxp3_counts.index) & (basemean_df_all[cols] > basemean_co).all(axis=1)
        v = wald_df_all[col].loc[idx].dropna()
        bins = pd.qcut(v, n, labels=False)
        bins.name = 'bin'
        
        # wald_df_by_bin = pd.concat([bins, foxp3_counts.loc[bins.index]], axis=1)
        v['bin'] = bins
        foxp3_count_by_bin = pd.concat([bins, foxp3_counts.loc[bins.index]], axis=1)
        xs = foxp3_count_by_bin.groupby('bin').mean(0).index
        ys = foxp3_count_by_bin.groupby('bin').mean(0).values
        xs -= np.median(xs)
        plt.plot(xs, ys, label=col_to_name[col],
                color = colors[c])
    plt.legend(bbox_to_anchor=(1, 1))
    plt.xlabel("Genes binned by Wald Statistics", labelpad=20)
    add_xaxis_labels('WT-fav.', 'Tir1-fav.', plt.gca(), fontsize=10)
    plt.ylabel(f"# {factor} peaks")
    plt.title(f"{factor} peaks by experiment")    



def make_motif_enrichment(foxp3_to_closest_tss, wald_clusters, wald_load_dict, foxp3_motif_df, wald_key = 'Bulk+sc Rest'):
    xs = []
    ys = []
    sums = []
    for u in sorted(wald_clusters[wald_key].dropna().unique()):
        genesoi = wald_clusters[wald_key].index[wald_clusters[wald_key]==u]
        atac_peaks = []
        for key, gene in foxp3_to_closest_tss.items():
            if gene[0] in genesoi:
                atac_peaks.append(key)
        atac_peaks = arr(atac_peaks).astype(int)
        xs.append(foxp3_motif_df.iloc[atac_peaks].mean())
        sums.append(foxp3_motif_df.iloc[atac_peaks].sum())
        ys.append(wald_load_dict[wald_key][u])
        
    rs = []
    motif_freq_df = pd.DataFrame(xs)
    for col in motif_freq_df:
        r, p = scipy.stats.pearsonr(motif_freq_df[col], ys)
        rs.append([col, r, p])

    rdf = pd.DataFrame(rs, columns = ['motif', 'r', 'p']).sort_values('r')
    return rdf, ys, motif_freq_df

def make_motif_counts(foxp3_to_closest_tss, wald_clusters, foxp3_motif_df, wald_key = 'Bulk+sc Rest'):
    # xs = []
    # ys = []
    sums = []
    for u in sorted(wald_clusters[wald_key].dropna().unique()):
        genesoi = wald_clusters[wald_key].index[wald_clusters[wald_key]==u]
        atac_peaks = []
        for key, gene in foxp3_to_closest_tss.items():
            if gene[0] in genesoi:
                atac_peaks.append(key)
        atac_peaks = arr(atac_peaks).astype(int)
        sums.append(foxp3_motif_df.iloc[atac_peaks])
    
    return sums
