#!/usr/bin/Rscript

#Read count scaling
#---
#Script to filter out low read count genes and scale the TPM read counts using z-scores.


#define global variables
DIR = "C:/Users/dstra/Nextcloud/Documents/Computers/Code/Github repos/bioinformatics analyses/Machine learning/Immune cell classification" #working directory
COUNT_CUTOFF = 10 #mean read count cutoff

#set working directory
setwd(DIR)

#load data
read_counts = read.csv("Data/TPM_gene_counts.csv", row.names = 1) #read count matrix


#Gene read filtering and scaling
#---
filtered_counts = read_counts[rowMeans(read_counts[]) > COUNT_CUTOFF,]

#scale the read counts using z-scores
scaled_counts = t(scale(t(filtered_counts), center = FALSE))

#save results as CSV
write.csv(scaled_counts, "Data/TPM_scaled_counts.csv", row.names = TRUE)
