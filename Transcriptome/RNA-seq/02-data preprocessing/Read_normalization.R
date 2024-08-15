#!/usr/bin/Rscript

#Read count normalization-DESeq2's median of ratios
#---
#Script to normalize the integer read counts for between-sample effects using DESeq2's median of ratios method.

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


#Read normalization-median of ratios (DESeq2)
#---
dds = DESeqDataSetFromMatrix(countData = read_counts, #create DESeq2 data object
	                         colData = metadata,
                             design = ~ 1)
dds = estimateSizeFactors(dds) #estimate size factors for each sample
norm_counts = counts(dds, normalized = TRUE) #create normalized read count matrix
norm_counts = data.frame(norm_counts) #convert to data frame

#save data frame as csv
write.csv(norm_counts, "norm_read_counts.csv", row.names = TRUE)
