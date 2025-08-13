import scanpy as sc
from copy import deepcopy

def normalize_adata(adata):
    adata.X = deepcopy(adata.layers['counts'].copy())
    sc.pp.normalize_total(adata)  
    sc.pp.log1p(adata)
    sc.pp.scale(adata, zero_center=True, max_value = 10)
    adata.layers['lognorm'] = adata.X.copy()
