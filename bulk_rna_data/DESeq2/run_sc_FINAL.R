source("./deg_funcs.R")


conditions <- c('rest', 'active')
for(condition in conditions){
    countData <- read.csv(glue('./sc_counts_FINAL/count_{condition}_df.csv'), row.names=1)
    colData <- read.csv('./sc_counts_FINAL/colData.csv', row.names=1)
    colData$HI = '1'
    if(dim(countData)[2] != dim(colData)[1]){
        print(glue("Skipping {condition}"))
        next
    }
    process_conditions(countData, colData, 'D0_tir1', 'D0_wt', glue('scFINAL_{condition}_'), thresh=0, design='sample')
    process_conditions(countData, colData, 'D3_tir1', 'D3_wt', glue('scFINAL_{condition}_'), thresh=0, design='sample')
    process_conditions(countData, colData, 'D7_tir1', 'D7_wt', glue('scFINAL_{condition}_'), thresh=0, design='sample')
}
