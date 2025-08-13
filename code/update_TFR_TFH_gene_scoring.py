import scanpy as sc
import pandas as pd
import numpy as np

def score_TFH(adata):
    TFH_data = pd.read_csv('/Genomics/argo/users/gdolsten/pritlab/wei_tir1_project/GSE124883_Tissue_Matrix.txt',
                    sep='\t', index_col = 0)
    tfh = TFH_data.loc[:, TFH_data.columns[TFH_data.columns.str.contains("LN Tfh")]].mean(axis=1)
    tfr = TFH_data.loc[:, TFH_data.columns[TFH_data.columns.str.contains("LN Tfr")]].mean(axis=1)
    
    treg = TFH_data.loc[:, TFH_data.columns[TFH_data.columns.str.contains("LN Treg")]].mean(axis=1)
    # (tfh / treg)[((treg>50) | (tfr>50))].sort_values().loc['Bcl6']
    
    sc.tl.score_genes(adata,
        (tfh / treg)[((treg>50) | (tfh>50))].sort_values().iloc[-200:].index,
                      score_name='Tfh > Treg'
                     )
    
    
    sc.tl.score_genes(adata,
        (tfr / treg)[((treg>50) | (tfr>50))].sort_values().iloc[-200:].index,
                      score_name='Tfr > Treg'
                     )
    
    
    sc.tl.score_genes(adata,
        (tfh / treg)[((treg>50) | (tfh>50))].sort_values().iloc[:200].index,
                      score_name='Tfh < Treg'
                     )
    
    
    sc.tl.score_genes(adata,
        (tfr / treg)[((treg>50) | (tfr>50))].sort_values().iloc[:200].index,
                      score_name='Tfr < Treg'
                     )
    
    
    sc.tl.score_genes(adata,
        (tfh / tfr)[((tfr>50) | (tfh>50))].sort_values().iloc[-200:].index,
                      score_name='Tfh > Tfr'
                     )
    
    sc.tl.score_genes(adata,
        (tfr / tfh)[((tfr>50) | (tfh>50))].sort_values().iloc[-200:].index,
                      score_name='Tfr > Tfh'
                     )
    