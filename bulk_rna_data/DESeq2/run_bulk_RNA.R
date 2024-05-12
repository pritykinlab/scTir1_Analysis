source("./deg_funcs.R")


#### BULK RNA ####
countData <- read.csv('count_df.csv', row.names=1)
colData <- read.csv('col_df.csv', row.names=1)

threshold <- 0
process_conditions(countData, colData, 'SpleenLN.TIR1.Activated', 'SpleenLN.WT.Activated', 'wei_',  thresh=threshold, design='Fullname', min_mean_reads=20)
process_conditions(countData, colData, 'Thymus.TIR1.Developing', 'Thymus.WT.Developing', 'wei_', thresh=threshold, design='Fullname', min_mean_reads=20)
process_conditions(countData, colData, 'SpleenLN.TIR1.Resting', 'SpleenLN.WT.Resting', 'wei_', thresh=threshold, design='Fullname', min_mean_reads=20)
process_conditions(countData, colData, 'SpleenLN.TIR1.Activated', 'SpleenLN.TIR1.Resting', 'wei_', thresh=threshold, design='Fullname')

#### BULK RNA Round 2 ####
countData <- read.csv('wei_bulk_round2_count_df.csv', row.names=1)
colData <- read.csv('wei_bulk_round2_col_df.csv', row.names=1)

countData <- countData[, colnames(countData) != 'tir1_d3_rep1']
colData <- colData[rownames(colData) != 'tir1_d3_rep1', ]

# colData$tmp = '1'
process_conditions(countData, colData, 'tir1_d7', 'wt_d7', 'weiNew_',  thresh=threshold, design='Celltype', min_mean_reads=20)
process_conditions(countData, colData, 'tir1_d3', 'wt_d3', 'weiNew_', thresh=threshold, design='Celltype', min_mean_reads=20)


#### Reintroduction Data ####
countData <- read.csv('reintroduction_count_df.csv', row.names=1)
colData <- read.csv('reintroduction_col_df.csv', row.names=1)
process_conditions(countData, colData, 'aKIKO.0.0.n', 'aKIKO.0.0.p',  'weiReintro_', thresh=threshold, design='Fullname', min_mean_reads=20)
process_conditions(countData, colData, 'aKIKO.3.0.n', 'aKIKO.3.0.p', 'weiReintro_', thresh=threshold, design='Fullname', min_mean_reads=20)
process_conditions(countData, colData, 'aKIKO.7.0.n', 'aKIKO.7.0.p',   'weiReintro_', thresh=threshold, design='Fullname', min_mean_reads=20)
process_conditions(countData, colData, 'aKIKO.14.0.n', 'aKIKO.14.0.p', 'weiReintro_', thresh=threshold, design='Fullname', min_mean_reads=20)
process_conditions(countData, colData, 'aKIKO.28.0.n', 'aKIKO.28.0.p',  'weiReintro_', thresh=threshold, design='Fullname', min_mean_reads=20)
process_conditions(countData, colData, 'nKIKO.0.0.n', 'nKIKO.0.0.p',  'weiReintro_', thresh=threshold, design='Fullname', min_mean_reads=20)
process_conditions(countData, colData, 'nKIKO.3.0.n', 'nKIKO.3.0.p',  'weiReintro_', thresh=threshold, design='Fullname', min_mean_reads=20)
process_conditions(countData, colData, 'nKIKO.7.0.n','nKIKO.7.0.p', 'weiReintro_', thresh=threshold, design='Fullname', min_mean_reads=20)
process_conditions(countData, colData, 'nKIKO.14.0.n', 'nKIKO.14.0.p',  'weiReintro_', thresh=threshold, design='Fullname', min_mean_reads=20)
process_conditions(countData, colData, 'nKIKO.28.0.n', 'nKIKO.28.0.p',  'weiReintro_', thresh=threshold, design='Fullname', min_mean_reads=20)


#### Joris Data ####
countData <- read.csv('joris_count_df.csv', row.names=1)
colData <- read.csv('joris_col_df.csv', row.names=1)
colData$tmp = '1'
process_conditions(countData, colData, 'activated_Foxp3_GFPKO', 'activated_Treg_GFP_DTR_WT', 'joris_', thresh=threshold, design='Celltype', min_mean_reads=20)
process_conditions(countData, colData, 'resting_Foxp3_GFPKO', 'resting_Treg_GFP_DTR_WT', 'joris_', thresh=threshold, design='Celltype', min_mean_reads=20)
process_conditions(countData, colData, 'CD73__Thymic_Foxp3_GFPKO', 'CD73__Thymic_Treg', 'joris_', thresh=threshold, design='Celltype', min_mean_reads=20)
