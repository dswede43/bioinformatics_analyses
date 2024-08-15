#!/usr/bin/Rscript

#Independent filtering (DESeq2)
#---
#This script completes independent filtering using DESeq2 to determine the normalized
#read count threshold for lowly expressed genes (genes that are unlikely to be
#differentially expressed). Genes below this threshold can be filtered out during
#differential expression analysis to reduce the multiple comparisons problem and
#improve power of statistical tests. Independent filtering is handled automatically
#by the DESeq2 package.
#(https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

#define global variables
read_count_files = c('read_count_matrix_file_name.csv','metadata_file_name.csv')
factors = c('categorical','variable','names')

#load packages
library(DESeq2)
source("https://raw.githubusercontent.com/dswede43/bioinformatics_analyses/main/Functions/load_read_count_data.R")


#Load the data
#---
read_count_data = load_read_count_data(read_count_files[1], read_count_files[2], factors)

#read count matrix
read_counts = read_count_data[[1]]

#experiment metadata
metadata = read_count_data[[2]]


#Independent filtering
#---
#define the generalized linear regression model
model = ~ your + independent + variables

#create the DESeq2 object
dds = DESeqDataSetFromMatrix(countData = read_counts,
							 colData = metadata,
							 design = model)

#complete differential expression analysis
dds = DESeq(dds)
DE_results = results(dds, alpha = 0.05)

#obtain normalized count threshold from independent filtering
metadata(DE_results)$filterThreshold

#plot a visual representation of independent filtering
plot(metadata(DE_results)$filterNumRej, 
    type = "b", ylab = "number of rejections",
    xlab = "quantiles of filter")
lines(metadata(DE_results)$lo.fit, col = "red")
abline(v = metadata(DE_results)$filterTheta)
