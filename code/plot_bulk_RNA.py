import sys
import numpy as np
from aux_functions import make_order_and_cluster_custom
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore
import seaborn as sns
from plotting_functions import init_subplots
from adjustText import adjust_text

import scipy
import scipy.stats

def format_pvalue(pval):
    if pval > 0.05:
        return 'NS'
    else:
        return f'{pval:.1e}'

import matplotlib.patheffects as path_effects
def plot_bulk_fc_fc_and_cdf_plot(comparisons, lfc_df_all, pval_df_all, basemean_df_all, PCO = .05, ylimdict = None,
                                 xlimdict = {}, 
                                        first_label = 'Tir1', second_label = 'GFPKO', basemean_co = 1,
                                        xlabel = 'TIR1 / WT', ylabel = 'KO GFP / WT GFP', lim=1, titledict = {},
                                        ):
    n = len(comparisons)
    fig = plt.figure(figsize=(4*n, 8))
    gs = gridspec.GridSpec(3, n, width_ratios=[1]*n, height_ratios=[1, .5, .5], wspace=.6, hspace=.6,)
    # The main heatmap is 9 times wider than the annotation heatmap
    for c, (cond, (comp_tir1, comp_ko)) in enumerate(comparisons.items()):
        z= np.linspace(-3, 3)
        if c == 0:
            ax = plt.subplot(gs[0, c])
            og_ax = ax
        else:
            ax = plt.subplot(gs[0, c], sharex=og_ax#, sharey=og_ax
            )
        high_tir1 = basemean_df_all[comp_tir1] > basemean_co
        high_ko = (basemean_df_all[comp_ko] > basemean_co)
        basemean_idx = high_tir1 | high_ko
        high_basemean_genes = basemean_df_all.index[basemean_idx]
        # return high_basemean_genes
        lfc_df, pval_df = lfc_df_all.loc[high_basemean_genes], pval_df_all.loc[high_basemean_genes]

        tir1, ko = lfc_df[comp_tir1], lfc_df[comp_ko]
        bad = np.isnan(tir1) | np.isnan(ko)
        r = scipy.stats.pearsonr(tir1[~bad], ko[~bad])[0].round(2)
        plt.text(1, 0, f'r = {r}', fontsize=12, ha='right', va='bottom', transform=ax.transAxes)
        pval_tir1, pval_ko = pval_df[comp_tir1], pval_df[comp_ko]
        idx_dict = {
            'neither' : (pval_tir1 > PCO) & (pval_ko > PCO),
            'tir1 up' : (tir1 > 0) & (pval_tir1 < PCO) & (pval_ko > PCO),
            'ko up' : (ko > 0) & (pval_tir1 > PCO) & (pval_ko < PCO),
            'tir1 down' : (tir1 < 0) & (pval_tir1 < PCO) & (pval_ko > PCO),
            'ko down' : (ko < 0) & (pval_tir1 > PCO) & (pval_ko < PCO),
            'up and down' : (tir1 > 0) & (ko < 0) & (pval_tir1 < PCO) & (pval_ko < PCO),
            'down and up' : (tir1 < 0) & (ko > 0) & (pval_tir1 < PCO) & (pval_ko < PCO),
            'both up' : (tir1 > 0) & (ko > 0) & (pval_tir1 < PCO) & (pval_ko < PCO),
            'both down' : (tir1 < 0) & (ko < 0) & (pval_tir1 < PCO) & (pval_ko < PCO),
        }
        color_dict = {
            'both up' : '#C5001D',
            'tir1 up' : '#FF6633',
            'ko up' : '#FF9999',
            'both down' : '#3366CC',
            'tir1 down' : '#66CCFF',
            'ko down' : '#9999FF',
            'neither' : '#EEEEEE',
            'up and down' : 'green',
            'down and up' : 'gold',
        }
        comp1, comp2 = cond.split("_vs_")

        for label, idx in idx_dict.items():
            color = color_dict[label]
            if 'ko' in label:
                color = '#EEEEEE'
            plt.scatter(tir1[idx], ko[idx], s=2, label=label, color=color, zorder = 3)
            plt.title(titledict.get(cond, cond))
            plt.xlabel(f"{comp1}")
            plt.ylabel(f"{comp2}")
            z = np.linspace(-5, 5)
            plt.plot(z, z, color='lightgray', linestyle='--', zorder=1)
            plt.ylim(ylimdict[cond])
            plt.yticks([ylimdict[cond][0], 0, ylimdict[cond][1]])
            plt.xlim(xlimdict.get(cond, [-2, 2]))
            plt.xticks([xlimdict.get(cond, [-2, 2])[0], 0, xlimdict.get(cond, [-2, 2])[1]])
        
        # Make CDF plot of comp using DEGs from other comp

        label_dict = {
            'tir1 up' : 'Up',
            'tir1 down' : 'Down',
            'ko up' : 'Up',
            'ko down' : 'Down',
            'neither' : 'NS',
        }
        comparison_dict = {
            f'{second_label} Subset' : [comp_ko, comp_tir1, comp2, comp1, 'ko', 'tir1'],
            f'{first_label} Subset' : [comp_tir1, comp_ko, comp1, comp2, 'tir1', 'ko'],
        }
        for _, (comp_dict_key, (comp_to_subset, comp_for_lfcs, subset_label, 
                                LFC_label, subset_label_for_labeldict, lfc_label_for_labeldict)) in enumerate(comparison_dict.items()):
            ax_comp = plt.subplot(gs[_ + 1, c])

            deg_dict = {
                subset_label_for_labeldict + ' down' : (pval_df[comp_to_subset] < PCO) & (lfc_df[comp_to_subset] < 0),
                subset_label_for_labeldict + ' up' : (pval_df[comp_to_subset] < PCO) & (lfc_df[comp_to_subset] > 0),
                'neither' : (pval_df[comp_to_subset] > PCO),
            }
            # if 'Activ' in cond:
                # print((pval_df[comp_to_subset] < PCO).sum())
                # for key, val in deg_dict.items():
                    # print(val.sum())
            baseline = lfc_df.loc[deg_dict['neither'], comp_for_lfcs]
            for c2, (label, degs) in enumerate(deg_dict.items()):
                vals = lfc_df.loc[degs, comp_for_lfcs]
                if len(vals) != 0:                
                    pval = scipy.stats.ks_2samp(vals, baseline)[1]
                else:
                    pval = 1
                color_label = label
                if label == 'neither':
                    sns.ecdfplot(vals, ax=ax_comp, label=f'{subset_label} DEGs {label_dict[label]};\nn={len(vals)}', 
                            color=color_dict[color_label], lw = 1.5,)
                else:
                    sns.ecdfplot(vals, ax=ax_comp, label=f'{subset_label} DEGs {label_dict[label]};\nn={len(vals)}; p={format_pvalue(pval)}', 
                            color=color_dict[color_label], lw = 1.5,)
                    
            plt.sca(ax_comp)
            plt.legend(bbox_to_anchor=(0, 1), loc='upper left', fontsize=6, frameon=False)
            ax_comp.set_xlabel(f"RNA LFC, {LFC_label}")

            ax_comp.set_ylabel("CDF")
            ax_comp.set_yticks([0, 0.5, 1])
            ax_comp.set_ylim([0, 1])
            ax_comp.set_title(f'{subset_label} DEGs')

            ax_comp.set_xlim([-lim, lim])
            ax_comp.set_xticks([-lim, 0, lim])
            if 'Tir1' in comp_dict_key:
                if cond == 'Thymic':
                    ax_comp.set_xlim([-lim, lim])
                    ax_comp.set_xticks([-lim, 0, lim])
                elif cond == 'Resting':
                    ax_comp.set_xlim([-lim, lim])
                    ax_comp.set_xticks([-lim, 0, lim])
                elif cond == 'Activated':
                    ax_comp.set_xlim([-lim, lim])
                    ax_comp.set_xticks([-lim, 0, lim])
            elif second_label in comp_dict_key:
                if cond == 'Thymic':
                    ax_comp.set_xlim([-lim, lim])
                    ax_comp.set_xticks([-lim, 0, lim])
                elif cond == 'Resting':
                    ax_comp.set_xlim([-lim, lim])
                    ax_comp.set_xticks([-lim, 0, lim])
                elif cond == 'Activated':
                    ax_comp.set_xlim([-lim, lim])
                    ax_comp.set_xticks([-lim, 0, lim])
    return fig


def plot_bulk_fc_fc_plot(comparisons, lfc_df_all, pval_df_all, basemean_df_all, PCO = .05, xlimdict = None, ylimdict = None, 
                                        basemean_co = 1, titledict = {}, label_x = 'tir1', label_y = 'ko',
                                        genes_to_label_dict = {}, labels_to_plot=None):
    n = len(comparisons)
    fig = plt.figure(figsize=(4*n, 8), dpi = 200)
    gs = gridspec.GridSpec(3, n, width_ratios=[1]*n, height_ratios=[1, .5, .5], wspace=.6, hspace=.6,)
    # The main heatmap is 9 times wider than the annotation heatmap
    for c, (cond, (comp_tir1, comp_ko)) in enumerate(comparisons.items()):
        z= np.linspace(-3, 3)
        if c == 0:
            ax = plt.subplot(gs[0, c])
        else:
            ax = plt.subplot(gs[0, c])
        plt.sca(ax)
        high_tir1 = basemean_df_all[comp_tir1] > basemean_co
        high_ko = (basemean_df_all[comp_ko] > basemean_co)
        basemean_idx = high_tir1 | high_ko
        high_basemean_genes = basemean_df_all.index[basemean_idx]
        # return high_basemean_genes
        lfc_df, pval_df = lfc_df_all.loc[high_basemean_genes], pval_df_all.loc[high_basemean_genes]

        tir1, ko = lfc_df[comp_tir1], lfc_df[comp_ko]
        bad = np.isnan(tir1) | np.isnan(ko)
        r = scipy.stats.pearsonr(tir1[~bad], ko[~bad])[0].round(2)
        plt.text(1, 0, f'r = {r}', fontsize=12, ha='right', va='bottom', transform=ax.transAxes)
        pval_tir1, pval_ko = pval_df[comp_tir1], pval_df[comp_ko]
        idx_dict = {
            'neither' : (pval_tir1 > PCO) & (pval_ko > PCO),
            f'{label_x} up' : (tir1 > 0) & (pval_tir1 < PCO) & (pval_ko > PCO),
            f'{label_y} up' : (ko > 0) & (pval_tir1 > PCO) & (pval_ko < PCO),
            f'{label_x} down' : (tir1 < 0) & (pval_tir1 < PCO) & (pval_ko > PCO),
            f'{label_y} down' : (ko < 0) & (pval_tir1 > PCO) & (pval_ko < PCO),
            'up and down' : (tir1 > 0) & (ko < 0) & (pval_tir1 < PCO) & (pval_ko < PCO),
            'down and up' : (tir1 < 0) & (ko > 0) & (pval_tir1 < PCO) & (pval_ko < PCO),
            'both up' : (tir1 > 0) & (ko > 0) & (pval_tir1 < PCO) & (pval_ko < PCO),
            'both down' : (tir1 < 0) & (ko < 0) & (pval_tir1 < PCO) & (pval_ko < PCO),
        }
        color_dict = {
            'both up' : '#C5001D',
            f'{label_x} up' : '#FF6633',
            f'{label_y} up' : '#FF9999',
            'both down' : '#3366CC',
            f'{label_x} down' : '#66CCFF',
            f'{label_y} down' : '#9999FF',
            'neither' : '#EEEEEE',
            'up and down' : 'green',
            'down and up' : 'gold',
        }
        zorder_dict = {
            'both up' : 1.4,
            f'{label_x} up' : 1.4,
            f'{label_y} up' : 1.3,
            'both down' : 1.4,
            f'{label_x} down' : 1.4,
            f'{label_y} down' : 1.3,
            'neither' : 1.3,
            'up and down' : 1.3,
            'down and up' : 1.3,
        }
        comp1, comp2 = cond.split("_vs_")
        if labels_to_plot is None:
            labels_to_plot = ['both up', 'both down', f'{label_x} up', f'{label_x} down'];
        label_count = 0
        for label, idx in idx_dict.items():
            color = color_dict[label]
            if 'ko' in label:
                color = '#EEEEEE'
            if label in labels_to_plot:
                plt.text(0, 1-.1*label_count, f'n = {idx.sum()}', fontsize=12, ha='left', va='top', transform=ax.transAxes, 
                         color=color)
                label_count += 1
            plt.scatter(tir1[idx], ko[idx], s=12, color=color, zorder = zorder_dict[label], edgecolor=None, linewidth=0)
            plt.scatter([], [], s=10, color=color, label=label)
            plt.title(titledict.get(cond, cond))
            plt.xlabel(f"{comp1}")
            plt.ylabel(f"{comp2}")
            z = np.linspace(-50, 50)
            plt.plot(z, z, color='lightgray', linestyle='--', zorder=1)
        plt.xlim(xlimdict[cond])
        plt.ylim(ylimdict[cond])
        plt.xticks([xlimdict[cond][0], 0, xlimdict[cond][1]])
        plt.yticks([ylimdict[cond][0], 0, ylimdict[cond][1]])
        genes_to_label = pval_tir1.index[idx_dict['both up'] | idx_dict['both down']]
        texts = []
        for gene in genes_to_label:
            if gene in high_basemean_genes:
                x, y = tir1.loc[gene], ko.loc[gene]
                if gene == 'Cxcr6':
                    print(x, y)
                t = plt.text(x, y, gene, fontsize=6, ha='right', va='bottom', color='black', zorder=20)
                plt.scatter(x, y, s=12, color='None', zorder=20, edgecolor='black', linewidth=.5)
                texts.append(t)
        adjust_text(texts, ax=plt.gca(), arrowprops=dict(arrowstyle="-", color='black', lw=0.15, alpha=.8), 
                    zorder=3.5, explode_radius = 150
        )

        if c == 0:
            plt.legend(bbox_to_anchor=(0, -.3), loc='upper left', fontsize=6, frameon=True, ncol=12)
    return fig







def make_fc_fc_plots(lfc_df, N_COLS_TAMOX, A_COLS_TAMOX, n_basecol = 'weiReintro_nKIKO.0.0.n_vs_nTreg.nan.nan_thresh=0',
            a_basecol = 'weiReintro_aKIKO.0.0.n_vs_aTreg.nan.nan_thresh=0'):
    
    N_COLS = len(N_COLS_TAMOX)
    z = np.linspace(-10, 10)
    fig, axs = init_subplots(len(N_COLS_TAMOX), 1)
    for c, col in enumerate(N_COLS_TAMOX):
        plt.sca(axs[c])
        name = '\n'.join(col.split("_vs_"))
        basename = '\n'.join(n_basecol.split("_vs_"))
        
        if col == n_basecol:
            continue
        xs = lfc_df.copy().loc[:, n_basecol]
        ys = lfc_df.copy().loc[:, col]
        plt.scatter(xs, ys, s = 2, zorder=4)
        plt.xlabel(name)
        plt.ylabel(basename)
        plt.title("Resting")
        plt.tight_layout()
        plt.xlim([-10, 10])
        plt.ylim([-10, 10])
        plt.plot(z, z, color='black')
        plt.grid(False)
        
    z = np.linspace(-10, 10)
    fig, axs = init_subplots(len(A_COLS_TAMOX), 1)
    for c, col in enumerate(A_COLS_TAMOX):
        plt.sca(axs[c])
        name = '\n'.join(col.split("_vs_"))
        basename = '\n'.join(a_basecol.split("_vs_"))
        
        if col == a_basecol:
            continue
        xs = lfc_df.copy().loc[:, a_basecol]
        ys = lfc_df.copy().loc[:, col]
        plt.scatter(xs, ys, s = 2, zorder=4)
        plt.xlabel(name)
        plt.ylabel(basename)
        plt.title("Active")
        plt.tight_layout()
        plt.xlim([-10, 10])
        plt.ylim([-10, 10])
        plt.plot(z, z, color='black')
        plt.grid(False)

def make_cdf_plots(lfc_df, N_COLS, ax, basecol = 'weiReintro_nKIKO.0.0.p_vs_nTreg.nan.nan_thresh=0', xlim=[-4, 4]):
    plt.sca(ax)
    colors = sns.color_palette("coolwarm", n_colors=len(N_COLS))
    for c, col in enumerate(N_COLS):
        name = '\n'.join(col.split("_vs_"))
        sns.ecdfplot(lfc_df.copy().loc[:, col], color=colors[c], label = name)
    plt.xlabel(name)
    plt.title("Resting")
    plt.tight_layout()
    plt.xlim(xlim)
    plt.legend(fontsize=6)

import matplotlib.gridspec as gridspec
import matplotlib.patches as patches

def plot_bulk_rna_heatmap_and_clusters(mat, masked_mat, lfc_df_all, pval_df_all, cluster, annot_cols, o, genes_to_plot, color_dict, color_quant_dict,
                                       cluster_max = 4, cbar_label=r'$ \mathrm{Log}_{2}$FC'):
    fig = plt.figure(figsize=(6, 6))
    gs = gridspec.GridSpec(1, 10, width_ratios=[1, 1, 1, 1, 1, 1, 1, 1, 1, 8],
                        )  # The main heatmap is 9 times wider than the annotation heatmap

    gs.update(wspace=0.2)

    # Plot the main heatmaps
    ax0 = plt.subplot(gs[-1])
    sns.heatmap(mat.iloc[o],
                cmap='bwr', vmin=-cluster_max, vmax=cluster_max,
                ax=ax0, cbar=False, yticklabels=False, xticklabels=True,)
    ax0.set_title(r'$ \mathrm{Log}_{2}$FC: Tir1 ÷ WT')

    # # Plot the annotation heatmaps

    axes = []
    previous_ax = None
    ax0_annotate = plt.subplot(gs[-2], sharey=ax0)

    cluster_df = pd.DataFrame(cluster, index=mat.index, columns=['Cluster'])
    # remap_color_dict = dict(zip(color_dict.keys(), np.arange(len(color_dict.keys()))))
    t = np.asarray([color_quant_dict[x[0]] for x in cluster_df.iloc[o].values])
    sns.heatmap(t[:, None], ax=ax0_annotate,
                cbar=False, cmap='coolwarm', vmin=-cluster_max, vmax=cluster_max, linewidth=0.005*0, linecolor='white')
    ax0_annotate.set_xticks([])
    ax0_annotate.set_xlabel("", rotation=90)

    previous_ax = ax0_annotate
    axes.append(ax0_annotate)

    degs_in_order = mat.iloc[o].index
    # for i, (key, col) in enumerate(annot_cols.items()):
    #     i = i+1
    #     ax0_annotate = plt.subplot(gs[i], sharey=ax0)
        
    #     # Remove spacing between this subplot and the previous one
    #     if previous_ax:
    #         box1 = previous_ax.get_position()
    #         box2 = ax0_annotate.get_position()
    #         ax0_annotate.set_position([box1.x1+.005, box2.y0, box2.width, box2.height])

    #     # Assuming annot_data1 contains your annotation data for the first heatmap
    #     annot_lfcs = np.sign(lfc_df_all.loc[degs_in_order, [col]]) * (pval_df_all.loc[degs_in_order, [col]] < .2)
    #     cbar = sns.heatmap(annot_lfcs, ax=ax0_annotate, 
    #                 cbar=False, cmap='bwr', vmin=-2, vmax=2, linewidth=0.8, linecolor='white')
        
    #     if i != 0:
    #         ax0_annotate.tick_params(labelleft=False, left=False) 
    #     else:
    #         tmp = ax0_annotate
            
    #     ax0_annotate.set_xticks([])
    #     if i == 0:
    #         ax0_annotate.set_xlabel("", rotation=90)
    #     else:
    #         ax0_annotate.set_xlabel(key.split(": ")[1], rotation=90)
        
    #     previous_ax = ax0_annotate
    #     axes.append(ax0_annotate)

    # for c, ax in enumerate(axes):
    #     box1 = ax.get_position()
    #     ax.set_position([box1.x0+.055, box1.y0, box1.width, box1.height])

    #     key = list(annot_cols.keys())[c-1]
    #     if 'Joris' in key:
    #         color = 'green'
    #     elif 'Saka' in key:
    #         color = 'red'
    #     else:
    #         color = 'black'

    #     if c > 0:
    #         rect = patches.Rectangle((box1.x0+.055, box1.y0 - box1.width-.01), 
    #                                 box1.width, box1.width, linewidth=1, 
    #                                 edgecolor='black', facecolor=color, transform=fig.transFigure)
    #         fig.add_artist(rect)
    #     else:
    #         plt.text(box1.x0+.055+box1.width*.9, box1.y0 - (box1.width + .01), "Dataset: ", transform=fig.transFigure, 
    #                     fontsize=12, ha='right', va='bottom')

    #     if ' r' in key:
    #         color = 'white'
    #     elif ' a' in key:
    #         color = 'gray'
    #     elif ' t' in key:
    #         color = 'black'

    #     if c > 0:
    #         rect = patches.Rectangle((box1.x0+.055, box1.y0 - 2*(box1.width + .01)), 
    #                                 box1.width, box1.width, linewidth=1, 
    #                                 edgecolor='black', facecolor=color, transform=fig.transFigure)
    #         fig.add_artist(rect)
    #     else:
    #         plt.text(box1.x0+.055+box1.width*.9, box1.y0 - 2*(box1.width + .01), "Status: ", transform=fig.transFigure, 
    #                     fontsize=12, ha='right', va='bottom')
    #     ax.xaxis.labelpad = 52

    #     plt.sca(ax)
    #     plt.grid(False)

    # for ax in [ax0]:
    #     ax.tick_params(labelleft=False, left=False) 
        
    #     plt.sca(ax)
    #     plt.grid(False)

    axes[0].set_yticks([])
    n = len(genes_to_plot)
    genes_to_plot = genes_to_plot[::-1]
    genes_to_plot = mat.iloc[o].index[np.isin(mat.iloc[o].index, genes_to_plot)]
    for c, i in enumerate(genes_to_plot):
        idx = np.where(mat.iloc[o].index == i)[0][0]
        width = 1/len(mat)/2
        location = 1 - (idx / len(mat) + width)
        ## Draw an arrow from text to location on heatmap:
        plt.sca(axes[0])
        textx, texty = -2, 1-(c+.5)/n
        cluster_val = cluster[o][idx]
        color_val = color_dict[cluster_val]
        plt.text(textx, texty, i, transform=axes[0].transAxes, ha='right', va='center', fontsize=8, color = 'black')
        plt.annotate('', xy=(textx, texty), xycoords='axes fraction', xytext=(-.2, location),
                    arrowprops=dict(arrowstyle="-", color = color_val, lw=.5))

    # # Make legend
    # legend_dict = {
    # 'green' : 'Joris',
    # 'red' : 'Sakaguchi',
    # 'black' : 'Active',
    # 'gray' : 'Resting',
    # 'white' : 'Thymic',
    # }
    # d = box1.width
    # shift = d*1.5
    # for c, (color, name) in enumerate(legend_dict.items()):
    #     rect = patches.Rectangle((1, .5 + shift*c), 
    #                             d, d, linewidth=1, 
    #                             edgecolor='black', facecolor=color, transform=fig.transFigure)
    #     plt.text(1+d*1.5, .5 + shift*c + d/2, name, transform=fig.transFigure, va='center')
    #     fig.add_artist(rect)

    cbar_ax = fig.add_axes([1, .2, 0.02, 0.45-.2])
    cbar = plt.colorbar(ax0.collections[0], cax=cbar_ax, )
    cbar_ax.set_ylabel(cbar_label)
    plt.grid(False)

    return fig, cluster_df
