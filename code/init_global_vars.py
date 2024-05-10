import matplotlib
import matplotlib.image
import numpy as np
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
import scanpy as sc
sc.set_figure_params(figsize=(5,5)) # no blurry figures allowed
sc.settings.verbosity = 4  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, color_map='viridis')

import sys
sys.path.append('../../code/')
sys.path.append('./code/')

from init_bulk_RNA import *
from init_bulk_ChIP import *
from init_bulk_ATAC import *
from aux_functions_scRNA import *
# from old_make_figures import *
from get_lfcs import generate_lfcs

PARSED_CHROMS = {'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX'}
PARSED_CHROMS_nochr = {
                '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                '11', '12', '13', '14', '15', '16', '17', '18', '19', 'X'}

import bbi 
import pybedtools as pbt
from aux_functions import *
