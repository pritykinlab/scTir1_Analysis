from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
from aux_functions_scRNA import aggregate_adata_by_obs_column
import pandas as pd
import numpy as np
import scipy.stats

def run_deseq(adata_dict, refill_padj = True, mean_cutoff = 5e-2):
    inference = DefaultInference(n_cpus=8)
    
    all_results_df = []
    for tissue in adata_dict:
        counts_df = aggregate_adata_by_obs_column(adata_dict[tissue], 'sample', 'counts', agg_func=np.sum).astype(int)
        mean_counts_df = aggregate_adata_by_obs_column(adata_dict[tissue], 'sample', 'counts', agg_func=np.mean)
        counts_df = counts_df[mean_counts_df.mean(axis=1) > mean_cutoff]
        print(f"{tissue} # genes : {counts_df.shape}")
        
        metadata = pd.DataFrame(index=counts_df.columns)
        metadata['genotype'] = metadata.index.str.split("_").str[1]
        
        dds = DeseqDataSet(
            counts=counts_df.T,  # transpose so rows = samples
            metadata=metadata,
            design = "~genotype",
            refit_cooks=True,
            inference=inference,
        )
        
        dds.deseq2()
    
        stat_res = DeseqStats(dds, contrast=('genotype', 'TIR1', 'WT'))
        stat_res.summary()
        results_df = stat_res.results_df
        results_df['tissue'] = tissue
        all_results_df.append(results_df)

    all_results_df = pd.concat(all_results_df, axis=0)
    deseq_l2fc_df = all_results_df.reset_index().pivot(columns='tissue', index='index', values='log2FoldChange')
    deseq_basemean_df = all_results_df.reset_index().pivot(columns='tissue', index='index', values='baseMean')
    if refill_padj:
        deseq_padj_df = all_results_df.reset_index().pivot(columns='tissue', index='index', values='pvalue')
        for col in deseq_padj_df:
            deseq_padj_df[col] =  scipy.stats.false_discovery_control(deseq_padj_df[col].fillna(1))
    else:
        deseq_padj_df = all_results_df.reset_index().pivot(columns='tissue', index='index', values='padj')
        
    return all_results_df, deseq_l2fc_df, deseq_padj_df, deseq_basemean_df






def run_deseq_matched(adata_dict, refill_padj = True, mean_cutoff = 5e-2):
    inference = DefaultInference(n_cpus=8)
    
    all_results_df = []
    all_dds = []
    for tissue in adata_dict:
        counts_df = aggregate_adata_by_obs_column(adata_dict[tissue], 'sample', 'counts', agg_func=np.sum).astype(int)
        mean_counts_df = aggregate_adata_by_obs_column(adata_dict[tissue], 'sample', 'counts', agg_func=np.mean)
        counts_df = counts_df[mean_counts_df.mean(axis=1) > mean_cutoff]
        
        metadata = pd.DataFrame(index=counts_df.columns)
        metadata['genotype'] = metadata.index.str.split("_").str[1]
        metadata['mouse'] = metadata.index.str.split("_").str[2]
        
        dds = DeseqDataSet(
            counts=counts_df.T,  # transpose so rows = samples
            metadata=metadata,
            design="~genotype + mouse",
            refit_cooks=True,
            inference=inference,
        )
        
        dds.deseq2()
    
        stat_res = DeseqStats(dds, contrast=('genotype', 'TIR1', 'WT'))
        stat_res.summary()
        results_df = stat_res.results_df
        results_df['tissue'] = tissue
        all_results_df.append(results_df)
        all_dds.append(dds)
    all_results_df = pd.concat(all_results_df, axis=0)
    deseq_l2fc_df = all_results_df.reset_index().pivot(columns='tissue', index='index', values='log2FoldChange')
    deseq_basemean_df = all_results_df.reset_index().pivot(columns='tissue', index='index', values='baseMean')
    if refill_padj:
        deseq_padj_df = all_results_df.reset_index().pivot(columns='tissue', index='index', values='pvalue')
        for col in deseq_padj_df:
            deseq_padj_df[col] =  scipy.stats.false_discovery_control(deseq_padj_df[col].fillna(1))
    else:
        deseq_padj_df = all_results_df.reset_index().pivot(columns='tissue', index='index', values='padj')
        
    return all_results_df, deseq_l2fc_df, deseq_padj_df, deseq_basemean_df, all_dds