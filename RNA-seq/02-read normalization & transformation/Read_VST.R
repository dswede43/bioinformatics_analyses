#!/usr/bin/Rscript

#Variance stabilizing transformation (VST)-DESeq2
#---
#Script to remove the heteroscedastic mean-variance relationship of read counts
#using the variance stabilizating transformation (VST) in DESeq2.

#define global variables
dir = "/path/to/directory"
read_count_file = 'your_read_counts.csv' #read count matrix file name
metadata_file = 'your_metadata.csv' #metadata file name
factors = c('your','categorical','variables') #factor variables

#set working directory
setwd(dir)

#load packages
library(dplyr)
library(DESeq2)
source("https://raw.githubusercontent.com/dswede43/bioinformatics_analyses/main/Functions/load_read_count_data.R")

#load data
data = load_read_count_data(read_count_file, metadata_file, factors)
read_counts = data[['read_counts']]
metadata = data[['metadata']]


#VST transformation (DESeq2)
#---
dds = DESeqDataSetFromMatrix(countData = read_counts, #create DESeq2 data object
	                         colData = metadata,
                             design = ~ 1)
vst_counts <- vst(dds, blind = FALSE) #obtain the vst transformed normalized read counts (blind = FALSE does not matter, model is ~1)
vst_counts <- assay(vst_counts)
vst_counts <- data.frame(vst_counts)

#save data frame as csv
write.csv(vst, "vst_read_counts.csv", row.names = TRUE)
