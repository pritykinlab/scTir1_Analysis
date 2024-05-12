library(DESeq2)
library(glue)


process_conditions <- function(countData, colData, c1, c2, prefix, thresh, design, min_mean_reads=20){
    indsoi <-  grepl(c1, colnames(countData)) | grepl(c2, colnames(countData))
    has_reads <- rowMeans(countData) > min_mean_reads
    countData <- countData[has_reads, ]
    design_var <- as.formula(paste0("~", design))

    countData <- countData[, indsoi]
    colData <- colData[indsoi, ]

    dds = DESeqDataSetFromMatrix(countData = countData, colData = colData, design=design_var)
    
    dds <- DESeq(dds)

    res = results(dds, contrast=c(design,  c1, c2), lfcThreshold = thresh )
    path <- glue('plots/{prefix}_dispersion_estimation_{c1}_vs_{c2}.png')
    png(path)
    plotDispEsts(dds)
    dev.off()

    path <- glue('plots/{prefix}_ma_plot_{c1}_vs_{c2}.png')
    png(path)
    plotMA(res, ylim=c(-2,2), xlim=c(.1, 1500))
    plotMA(res, )
    dev.off()
    path <- glue("output/{prefix}{c1}_vs_{c2}_thresh={thresh}.csv")
    write.table(res, path)
    sizefactors <- sizeFactors(dds)
    write.table(sizefactors, glue("sizefactors/{prefix}{c1}_vs_{c2}.csv"))
    return(res)
}

