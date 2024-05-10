import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
import numpy as np
import scipy
import scipy.stats

def plot_cdf_plots(columns_to_plot_dict, peak_metadata_df, full_lfc_df, baseline='NS'):
    statuses = {'cell_is_resting', 'cell_is_active'}
    for name, factors in columns_to_plot_dict.items():
        vdict = {}
        fig, axs = plt.subplots(2, 3, figsize=(12, 6))
        axs = np.ravel(axs)
        for factor in factors:
            if factor not in peak_metadata_df.columns:
                raise ValueError(f"{factor} not in peak_metadata_df")
            peaks_with_factor = peak_metadata_df.index[peak_metadata_df[factor] > 0]
            c = 0
            for status in statuses:
                for _, u in enumerate(sorted(full_lfc_df['comparison'].unique())):
                    ax = axs[c]
                    comp_inds = (full_lfc_df['comparison']==u)
                    stat_inds = (full_lfc_df['active_status']==status)
                    indsoi = (comp_inds & stat_inds)
                    subdf = full_lfc_df[indsoi].copy()
                    peak_inds = (subdf['gene'].isin(peaks_with_factor))
                    with_factor = subdf.loc[peak_inds]['log2FC']
                    n = len(with_factor)
                    if baseline in factor:
                        color = 'black'
                        text = baseline
                        text = text + f'\n n={n};'
                        vdict[baseline] = with_factor
                    else:
                        p = scipy.stats.ks_2samp(vdict[baseline], with_factor)[1]
                        text = '_'.join(factor.split("_")[:2]) + "-up"
                        text = text + f'\n n={n}; p={p:.0e}'
                        color = None
                    sns.ecdfplot(with_factor, ax=ax, label=f'{text}', color=color)
                    ax.set_xlim([-.8, .8])
                    ax.legend(fontsize=8)
                    ax.set_title(f'{status}_{u}')
                    ax.set_xlabel("LFC, scATAC")
                    ax.set_ylabel("Fraction")
                    c+=1
        plt.tight_layout()
        fig.suptitle(name, va='bottom')



import pybedtools as pbt
import os
import sys
sys.path.append('../../code')
from motif_analysis import fimo_to_bed
import pandas as pd
from aux_functions import add_chr_to_bedtool, get_col
import os

def collect_from_fimo(motif_path, atac_peak_names, all_atac_peaks, motif_folder):
    # print("Collecting")
    motif_path_name = os.path.basename(motif_path)
    _, motif_identity, motif_name = motif_path_name.lower().split("_")    
    motif_col_name = f'{motif_name}_{motif_identity}'

    fimo_tsv_path = f'./{motif_path}/fimo.tsv'
    fimo_bedfile_path = f'motif_bedfiles/{motif_folder}/{motif_path_name}.bed'
    os.makedirs(f'motif_bedfiles/{motif_folder}', exist_ok=True)
    try:
        fimo_to_bed(fimo_tsv_path, fimo_bedfile_path)
    except:
        motif_count_df = pd.DataFrame()  
        motif_count_df.index = atac_peak_names
        motif_count_df[motif_col_name] = 0   
        return motif_count_df
    motif_counts = pbt.BedTool(fimo_bedfile_path)
    counts = get_col(all_atac_peaks.intersect(motif_counts, c=True), -1).astype(float)
    motif_count_df = pd.DataFrame()
    motif_count_df[motif_col_name] = counts
    motif_count_df.index = atac_peak_names
    # print("Done Collecting")
    
    return motif_count_df

import pandas as pd
from collections import defaultdict
from aux_functions import remove_chr_bedtool
def atac_peak_to_closest_gene(tss_df, atac_bedtool, k=1, maxdist=np.inf, drop_bad=True):
    tss_bedtool = pbt.BedTool.from_dataframe(tss_df.reset_index()[['chrom', 'start', 'end', 'gene_name']])
    atac_to_closest_genes = defaultdict(list)
    atac_bedtool = pbt.BedTool.from_dataframe(atac_bedtool.to_dataframe().iloc[:, :4])
    z = atac_bedtool.sort().closest(tss_bedtool.sort(), k=k, d=True, t='first')
    for row in z:
        d = row[-1]
        if d == '.':
            genename = '.'
        if int(d) > maxdist:
            genename=  '.'
        genename = row[7]
        if (genename == '.') and (drop_bad):
            continue
        atac_peak = row[3]
        atac_to_closest_genes[atac_peak].append(genename)
    for i, v in atac_to_closest_genes.items():
        if len(atac_to_closest_genes[i])==0:
            atac_to_closest_genes[i].append('.')
    return atac_to_closest_genes


def atac_peak_to_closest_gene_NEW(tss_df, atac_bedtool, k=1, maxdist=np.inf):
    tss_bedtool = pbt.BedTool.from_dataframe(tss_df.reset_index()[['chrom', 'start', 'end', 'gene_name']])
    atac_to_closest_genes = defaultdict(list)
    atac_bedtool = pbt.BedTool.from_dataframe(atac_bedtool.to_dataframe().iloc[:, :4])
    z = atac_bedtool.sort().closest(tss_bedtool.sort(), k=k, d=True, t='first')
    for row in z:
        d = row[-1]
        if d == '.':
            continue
        if int(d) > maxdist:
            continue
        genename = row[6]
        atac_peak = '_'.join(row[:3])
        atac_to_closest_genes[atac_peak].append(genename)
    for i, v in atac_to_closest_genes.items():
        if len(atac_to_closest_genes[i])==0:
            atac_to_closest_genes[i].append('.')
    return atac_to_closest_genes



