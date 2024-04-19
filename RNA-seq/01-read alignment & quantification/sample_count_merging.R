#!/usr/bin/Rscript

#Read count sample merging
#---
#This script merges the read quantifications from each sample into a single read count matrix.

#load packages
library(dplyr)

#define global variables
dir = "/media/sf_VM_share/RNA-seq_data/KneeOA/outputs/quants" #working directory

#set working directory
setwd(dir)


#Load the quantification files
#---
#define the sample names
samples = list.dirs(full.names = FALSE)
samples = samples[samples != ""]

#load in each quantification file
sample_counts = list()
for (sample in samples){
    sample_counts[[sample]] = read.table(paste0(sample, "/sample", sample, ".txt"), header = TRUE, sep = "\t")
}


#Merge read count data into matrix
#---
count_matrix = do.call(cbind, sample_counts)
count_matrix = dplyr::select(count_matrix, contains("X"))
colnames(count_matrix) = samples
rownames(count_matrix) = sample_counts[[1]]$Geneid

#save the matrix as CSV file
write.csv(count_matrix, paste0(dir, "/read_count_matrix.csv"), row.names = TRUE)
