import glob
import pybedtools as pbt
# from init_bulk_ChIP import all_peak_dict
import sys
sys.path.append('../../code/')
from aux_functions import get_col, remove_chr_bedtool
atac_dict = {}
for file in glob.glob('../../atac/processed/*thresh=0*.csv'):
    cond1, cond2 = file.split("/")[-1].split("_thresh")[0].split("_vs_")
    thresh = file.split(".csv")[0].split("thresh=")[1]
    atac_dict[f'{cond1}_vs_{cond2}_thresh={thresh}'] = pbt.BedTool(file)

atac_dict_up_and_down = {}
for file in glob.glob('../../atac/processed/*thresh=*.csv'):
    cond1, cond2 = file.split("/")[-1].split("_thresh")[0].split("_vs_")
    thresh = file.split(".csv")[0].split("thresh=")[1]
    bedup = pbt.BedTool(file).filter(lambda x: (float(x[4]) > 0) & (float(x[5]) < .05)).saveas()
    beddown = pbt.BedTool(file).filter(lambda x: (float(x[4]) < 0) & (float(x[5]) < .05)).saveas()
    bedns = pbt.BedTool(file).filter(lambda x: (float(x[5]) > .05)).saveas()
    atac_dict_up_and_down[f'{cond1}_{cond1}_vs_{cond2}_thresh={thresh}'] = bedup
    atac_dict_up_and_down[f'{cond2}_{cond1}_vs_{cond2}_thresh={thresh}'] = beddown
    atac_dict_up_and_down[f'NS_{cond1}_vs_{cond2}_thresh={thresh}'] = bedns

def get_index_from_bedtool_df(df):
    ind = df['chrom'].astype(str) + "_" + df['start'].astype(str) + "_" + df['end'].astype(str) 
    return ind



import pandas as pd
atac_peak_bedtool = atac_dict['Treg_actv_vs_Tcon_actv_thresh=0']
atac_peak_df = atac_peak_bedtool.to_dataframe()
ind = get_index_from_bedtool_df(atac_peak_df)
atac_peak_df.index = ind

# cols = list(atac_dict_up_and_down.keys())+list(all_peak_dict.keys())
# atac_metadata_df = pd.DataFrame(index=ind, columns=cols)
# for col in cols:
# 	atac_metadata_df[col] = 0

# for name, bedtool in atac_dict_up_and_down.items():
#     inds = get_index_from_bedtool_df(bedtool.to_dataframe())
#     atac_metadata_df.loc[inds, name] = 1
    
# for name, bedtool in all_peak_dict.items():
#     atac_in_bedtool = atac_peak_bedtool.intersect(remove_chr_bedtool(bedtool), u=True)
#     inds = get_index_from_bedtool_df(atac_in_bedtool.to_dataframe())
#     atac_metadata_df.loc[inds, name] = 1