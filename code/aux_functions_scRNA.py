import numpy as np
import pandas as pd
import scipy.stats
from scipy.stats import mannwhitneyu, gmean

import statsmodels
import statsmodels.stats
import statsmodels.stats.multitest

arr = np.asarray
def make_df_aggregated_by_key(adata, genesoi, key='sample', leiden='leiden', layer='', useLog=False, pc=0, agg_func=np.mean):
    rows = []
    for u in np.unique(adata.obs[leiden]):
        cluster_inds = adata.obs[leiden]==u
        subadata = adata[cluster_inds, genesoi]
        for w in np.unique(subadata.obs[key]):
            sample_inds = subadata.obs[key] == w
            val = agg_func(subadata[sample_inds].layers[layer], axis=0).copy()
            if type(val) != np.ndarray:
                val = np.ravel(arr(val))
            val = val + pc
            if useLog:
                val = np.log2(val)
            for gene, v in zip(genesoi, val):
                row = [v, gene, u, w, w]
                rows.append(row)
    wei_df = pd.DataFrame(rows, columns=['value', 'gene', leiden, 'sample', 'parsed_sample'])
    wei_df['parsed_sample'] = wei_df['parsed_sample'].str.replace("_rep[0-9]", '', regex=True)
    return wei_df



from scipy.stats import mannwhitneyu, gmean

## DESeq2-like size factor
def calculate_size_factor_for_each_cell(adata, raw_layer):
    counts = np.ravel(adata.layers[raw_layer].sum(axis = 1).squeeze())
    g_mean = gmean(counts)
    size_factor = counts/g_mean
    return size_factor

import random
import math
import numpy as np
import matplotlib.pyplot as plt
def init_subplots(n, rows, fgsz=(4, 4   )):
    cols = math.ceil(n/rows)
    fig, axs = plt.subplots(rows, cols, figsize=(cols*fgsz[0], rows*fgsz[1]))
    axs = np.ravel(axs)
    return fig, axs

def aggregate_obs_adata_by_obs_column(adata, num_key, key='sample'):
    res_df = pd.DataFrame()
    mus = []
    us = sorted(np.unique(adata.obs[key]))
    for u in us:
        indsoi = adata.obs[key] == u
        mu = adata[indsoi].obs[num_key].mean()
        mus.append(mu)
    res_df[num_key] = mus
    res_df.index = us
    return res_df


import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF

def compute_adata_obs_pvalues(adata, control_conditions, test_conditions, numerical_column, condition_column):
    baseline_idx = adata.obs[condition_column].isin(control_conditions)
    nonbaseline_idx = adata.obs[condition_column].isin(test_conditions)

    control_values = adata.obs.loc[baseline_idx, numerical_column]
    test_values = adata.obs.loc[nonbaseline_idx, numerical_column]
    ecdf = ECDF(control_values)
    ecdf_at_new_points = ecdf(test_values)
    return ecdf_at_new_points, control_values




def aggregate_obs_adata_by_obs_column_and_leiden(adata, num_key, key='sample', leiden_key='leiden'):
    results_df = {}
    for u in sorted(np.unique(adata.obs[leiden_key])):
        indsoi = adata.obs[leiden_key] == u
        subadata = adata[indsoi]
        res_df = aggregate_obs_adata_by_obs_column(subadata, num_key, key=key)
        results_df[u] = res_df
    return results_df

def aggregate_adata_by_obs_column(adata, key='sample', layer=''):
    res_df = pd.DataFrame()
    us = np.unique(adata.obs[key])
    for u in us:
        indsoi = adata.obs[key] == u
        mu = adata.layers[layer][indsoi].mean(axis=0)
        res_df[u] = np.ravel(mu)
        
        res_df.index = adata.var.index
    return res_df


def make_df_aggregated_by_key_with_sizefactors(adata, genesoi, key='sample', leiden='leiden', layer='', 
                                                useLog=False, pc=0, agg_func=np.mean, seed=0, max_cells=np.inf):
    random.seed(seed)
    rows = []
    for u in np.unique(adata.obs[leiden]):
        cluster_inds = adata.obs[leiden]==u
        subadata = adata[cluster_inds, genesoi]

        for w in np.unique(subadata.obs[key]):
            sample_inds = subadata.obs[key] == w
            subsubadata = subadata[sample_inds]
            n = len(subsubadata)
            if len(subsubadata) > max_cells:
                random_inds = random.sample(range(n), max_cells)
                subsubadata = subsubadata[random_inds]
            data_for_sample = subsubadata

            size_factors = calculate_size_factor_for_each_cell(data_for_sample, layer)
            scaled_data = data_for_sample.layers[layer]/(size_factors[:, None])
            val = agg_func(scaled_data, axis=0).copy()

            if type(val) != np.ndarray:
                val = np.ravel(arr(val))
            val = val + pc
            if useLog:
                val = np.log2(val)
            for gene, v in zip(genesoi, val):
                row = [v, gene, u, w, w]
                rows.append(row)
    wei_df = pd.DataFrame(rows, columns=['value', 'gene', leiden, 'sample', 'parsed_sample'])
    wei_df['parsed_sample'] = wei_df['parsed_sample'].str.replace("_rep[0-9]", '', regex=True)
    return wei_df





def DEG_test_using_pseudobulk_pearson_residual(adata, pearson_layer, raw_layer, frac=.1):
    pairwise_test_df = pd.DataFrame()
    for key in ['cell_is_resting', 'cell_is_active']:
        genes_with_high_counts_indicator = np.ravel(arr(adata.layers[raw_layer].sum(axis=0)) > adata.shape[0] * frac)
        genes_with_high_counts = adata.var.index[genes_with_high_counts_indicator]

        wei_df = make_df_aggregated_by_key(adata, genes_with_high_counts, key='sample', leiden=key, layer=pearson_layer)
        wei_df = wei_df[wei_df[key] == 1]
        wei_df['key'] = key
        genes = wei_df['gene'].unique()
        for comp1, comp2 in [['D3_tir1', 'D3_wt'], ['D7_tir1', 'D7_wt'], ['D7_tir1', 'D3_tir1'], ['D0_tir1', 'D0_wt']]:
            pairwise_tests = []
            for gene in genes:
                subdf = wei_df[wei_df['gene'] == gene].copy()
                tmp_df = subdf.set_index('parsed_sample')
                v1, v2 = tmp_df.loc[comp1, 'value'], tmp_df.loc[comp2, 'value']
                stat, pval = scipy.stats.ttest_ind(v1, v2)
                name = f'{comp1}_vs_{comp2}'
                pairwise_tests.append([name, stat, pval, gene, key])
            tmpdf = pd.DataFrame(pairwise_tests, columns = ['comparison', 'statistic', 'pval', 'gene', 'key'])
            tmpdf['padj'] = statsmodels.stats.multitest.fdrcorrection(tmpdf['pval'])[1]
            pairwise_test_df = pd.concat([pairwise_test_df, tmpdf], axis=0)
    return pairwise_test_df

def make_matrix_from_DEG_output_df_cluster_specific(pearson_stat, name_subset, value_key='padj', key_order=None):
    df = pd.DataFrame()
    if key_order is None:
        us = np.unique(pearson_stat['key'])
    else:
        us = key_order
    for u in us:
        rest_df = pearson_stat[pearson_stat['key'] == u]
        pearson_rest_df = rest_df[rest_df['gene'].isin(name_subset)].pivot(index='gene', columns='comparison', values=value_key)
        df = pd.concat([df, pearson_rest_df.add_suffix(f" Cluster {u}")], axis=1)
    return df

def make_matrix_from_DEG_output_df(pearson_stat, name_subset, value_key='padj'):
    rest_df = pearson_stat[pearson_stat['key'] == 'cell_is_resting']
    pearson_rest_df = rest_df[rest_df['gene'].isin(name_subset)].pivot(index='gene', columns='comparison', values=value_key)
    print((pearson_rest_df < .05).sum())

    active_df = pearson_stat[pearson_stat['key'] == 'cell_is_active']
    pearson_active_df = active_df[active_df['gene'].isin(name_subset)].pivot(index='gene', columns='comparison', values=value_key)
    return pd.concat([pearson_rest_df.add_suffix(" Rest"), pearson_active_df.add_suffix(" Actv")], axis=1)

import scipy
import scipy.stats
import scipy.stats.mstats
from scipy.stats.mstats import gmean
def DEG_test_using_pseudobulk_lfcs(adata, raw_layer, recalculate=True, frac = .1, genes_to_test=None):
    def cellwise_size_factors(adata):
            counts_per_cell = arr(adata.layers[raw_layer].sum(axis = 0)).squeeze()
            gmean_vals = gmean(counts_per_cell, axis=0)
            size_factors = np.ravel(counts_per_cell/gmean_vals)
            return size_factors

    if ('raw_counts_by_size_factors' not in adata.layers.keys()) or (recalculate==True):
        adata.layers['raw_counts_by_size_factors'] = arr(adata.layers[raw_layer].todense()) * cellwise_size_factors(adata)
    else:
        print("Using existing size factors & raw counts matrix")
    pairwise_test_df = pd.DataFrame()
    for key in ['cell_is_resting', 'cell_is_active']:
        pairwise_tests = []
        # genes_with_high_counts_indicator = np.ravel(arr(adata.layers[raw_layer].sum(axis=0)) > adata.shape[0] * frac)
        # genes_with_high_counts = adata.var.index[genes_with_high_counts_indicator]
        if genes_to_test is None:
            genes_to_test = adata.var.index
        wei_df = make_df_aggregated_by_key(adata,  genes_to_test, key='sample', leiden=key, 
                                            useLog=True, layer='raw_counts_by_size_factors', pc=1)
        wei_df = wei_df[wei_df[key] == 1]
        wei_df['key'] = key
        genes = wei_df['gene'].unique()
        for gene in genes:
            subdf = wei_df[wei_df['gene'] == gene].copy()
            for group in subdf.groupby(key):
                tmp_df = group[1].set_index('parsed_sample')
                for comp1, comp2 in [['D3_tir1', 'D3_wt'], ['D7_tir1', 'D7_wt'], ['D7_tir1', 'D3_tir1'], ['D0_tir1', 'D0_wt']]:
                    v1, v2 = tmp_df.loc[comp1, 'value'], tmp_df.loc[comp2, 'value']
                    stat, pval = scipy.stats.ttest_ind(v1, v2)
                    name = f'{comp1}_vs_{comp2}'
                    pairwise_tests.append([name, stat, pval, gene, key])
        tmpdf = pd.DataFrame(pairwise_tests, columns = ['comparison', 'statistic', 'pval', 'gene', 'key'])
        tmpdf['pval'] = tmpdf['pval'].fillna(1)
        tmpdf['padj'] = statsmodels.stats.multitest.fdrcorrection(tmpdf['pval'])[1]
        pairwise_test_df = pd.concat([pairwise_test_df, tmpdf], axis=0)
    return pairwise_test_df


def DEG_test_using_pseudobulk_pearson_residual_other_comparison(adata, pearson_layer, raw_layer, frac=.1):
    pairwise_test_df = pd.DataFrame()
    for key in ['cell_is_resting', 'cell_is_active']:
        pairwise_tests = []
        genes_with_high_counts_indicator = np.ravel(arr(adata.layers[raw_layer].sum(axis=0)) > adata.shape[0] * frac)
        genes_with_high_counts = adata.var.index[genes_with_high_counts_indicator]
        wei_df = make_df_aggregated_by_key(adata, genes_with_high_counts, key='sample', leiden=key, layer=pearson_layer)
        wei_df = wei_df[wei_df[key] == 1]
        wei_df['key'] = key
        genes = wei_df['gene'].unique()
        for gene in genes:
            subdf = wei_df[wei_df['gene'] == gene].copy()
            for group in subdf.groupby(key):
                for control in ['D3_wt', 'D7_wt', 'D0_wt', 'D0_tir1']:
                    group[1]['parsed_sample'].replace(control, 'Controls', inplace=True)
                tmp_df = group[1].set_index('parsed_sample')

                for comp1, comp2 in [['D3_tir1', 'Controls'], ['D7_tir1', 'Controls']]:
                    v1, v2 = tmp_df.loc[comp1, 'value'], tmp_df.loc[comp2, 'value']
                    stat, pval = scipy.stats.ttest_ind(v1, v2)
                    name = f'{comp1}_vs_{comp2}'
                    pairwise_tests.append([name, stat, pval, gene, key])
        tmpdf = pd.DataFrame(pairwise_tests, columns = ['comparison', 'statistic', 'pval', 'gene', 'key'])
        tmpdf['padj'] = statsmodels.stats.multitest.fdrcorrection(tmpdf['pval'])[1]
        pairwise_test_df = pd.concat([pairwise_test_df, tmpdf], axis=0)
    return pairwise_test_df


def compute_cutoff(v, n_sigmas):
    mu = np.nanmean(v)
    sigma = np.nanstd(v)
    cutoff = mu + (n_sigmas*sigma)
    return cutoff

def cutoff_vector(v, n_sigmas):
    cutoff = compute_cutoff(v, n_sigmas)
    cutoff_vector = v > cutoff
    return cutoff_vector

def make_obs_aggregated_by_key(adata, obssoi, key='sample', agg_func = np.mean, leiden_key = 'leiden'):
    rows = []
    for u in np.unique(adata.obs[leiden_key]):
        cluster_inds = adata.obs[leiden_key]==u
        subadata = adata[cluster_inds]
        for w in np.unique(subadata.obs[key]):
            sample_inds = subadata.obs[key] == w
            val = agg_func(subadata[sample_inds].obs[obssoi], axis=0)
            for gene, v in zip(obssoi, val):
                row = [v, gene, u, w, w.replace("rep*", "")]
                rows.append(row)
    wei_df = pd.DataFrame(rows, columns=['value', 'gene', 'cluster', 'sample', 'parsed_sample'])
    return wei_df


def compute_summary_statistics_of_adata(adata, key = 'sample'):
    fract_wt_dict = {}
    n_cell_dict = {}
    fract_d7_tir1_dict = {}
    fract_d3_tir1_dict = {}
    fract_d0_tir1_dict = {}

    for u in adata.obs['leiden'].unique():
        frac_wt = adata.obs[adata.obs['leiden'] == u]['is_wt'].mean()
        n_cells = (adata.obs['leiden'] == u).sum()/adata.obs['leiden'].value_counts().max()
        fract_wt_dict[u] = frac_wt
        fract_d7_tir1_dict[u] = adata.obs[adata.obs['leiden'] == u][key].str.contains("D7_tir1").mean()
        fract_d3_tir1_dict[u] = adata.obs[adata.obs['leiden'] == u][key].str.contains("D3_tir1").mean()

        fract_d0_tir1_dict[u] = adata.obs[adata.obs['leiden'] == u][key].str.contains("D0_tir1").mean()
        n_cell_dict[u] = n_cells
    adata.obs['frac_wt_in_cluster'] = [fract_wt_dict[x] for x in adata.obs['leiden']]
    adata.obs['frac_d7_tir1_in_cluster'] = [fract_d7_tir1_dict[x] for x in adata.obs['leiden']]
    adata.obs['frac_d3_tir1_in_cluster'] = [fract_d3_tir1_dict[x] for x in adata.obs['leiden']]

    adata.obs['frac_d0_tir1_in_cluster'] = [fract_d0_tir1_dict[x] for x in adata.obs['leiden']]
    adata.obs['n_cells_in_cluster'] = [n_cell_dict[x] for x in adata.obs['leiden']]    


