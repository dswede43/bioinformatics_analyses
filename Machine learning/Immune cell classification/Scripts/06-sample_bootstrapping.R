#!/usr/bin/Rscript

#Sample bootstrapping
#---
#Script to bootstrap RNA-seq samples to increase sample size and split into
#training and testing datasets for machine learning.


#Table of contents:
#---
#1. Bootstrap sample
#2. Train and test data splits


#load packages
library(dplyr)
library(tidyr)
library(reshape2)
library(jsonlite)

#define global variables
DIR = "C:/Users/dstra/Nextcloud/Documents/Computers/Code/Github repos/bioinformatics analyses/Machine learning/Immune cell classification" #working directory
BOOT_TRIALS = 500 #number of bootrapping trials
TRAIN_SIZE = 0.8 #roportion of data for training set

#set working directory
setwd(DIR)

#load data
read_counts = read.csv("Data/TPM_gene_counts.csv") #scaled read count matrix
hub_genes = fromJSON("Data/module_hub_genes.json") #WGCNA modules
metadata = read.csv("Data/metadata.csv") #study metadata
metadata = metadata[,c('sample','cell_group')]


#1. Bootstrap sample
#---
#define an empty array
genes = c()
for (module_name in names(hub_genes)){
	#append the genes from each module into the array
	genes = c(genes, hub_genes[[module_name]]$ensembl)
}

#filter the scaled read counts for gene modules
read_counts = filter(read_counts, ensembl %in% genes)

#transform the scaled read count data into long format
count_data = melt(read_counts, id = 'ensembl')
colnames(count_data) = c("ensembl","sample","read_counts")

#merge the count data with metadata
count_data = merge(count_data, metadata, by = 'sample')
count_data$cell_group = ifelse(count_data$cell_group == 'T_cell', 'T_cell', 'non_T_cell')

#define the levels from the group variable
groups = unique(count_data$cell_group)

#define an empty data frame
boot_data = data.frame()
for (gene in genes){
	#filter the data for the current gene
	gene_data = filter(count_data, ensembl == gene)

	for (group in groups){
		#filter the data for the current group level
		gene_group_data = filter(gene_data, cell_group == group)

		#define the bootstrap sample size
		sample_size = nrow(gene_group_data)
		ntrials = BOOT_TRIALS %/% sample_size
		sample_size = sample_size * ntrials

		#bootstrap the data (sample with replacement)
		boot_samples = sample(gene_group_data$read_counts, sample_size, replace = TRUE)

		#organize the bootstrap data
		gene_boot_data = data.frame(sample = "sample", ensembl = gene, read_counts = boot_samples, cell_group = group)
		gene_boot_data = rbind(gene_group_data, gene_boot_data)
		gene_boot_data$sample = paste0(group, "_", 1:nrow(gene_boot_data))

		#store the results
		boot_data = rbind(boot_data, gene_boot_data)
	}
}

#transform the bootstrap data into matrix format
boot_counts = pivot_wider(boot_data, names_from = ensembl, values_from = read_counts)


#2. Train and test data splits
#---
#determine the training dataset size
train_size = round(nrow(boot_counts) * TRAIN_SIZE, 0)

#set the seed
set.seed(42)

#determine the training dataset indices
train_index = as.numeric(sample(rownames(boot_counts), train_size, replace = FALSE))

#split the data into training and testing sets
train_set = boot_counts[train_index, ]
test_set = boot_counts[-train_index, ]

#save the results as CSV
write.csv(train_set, "Data/ML_train_data.csv", row.names = FALSE)
write.csv(test_set, "Data/ML_test_data.csv", row.names = FALSE)
