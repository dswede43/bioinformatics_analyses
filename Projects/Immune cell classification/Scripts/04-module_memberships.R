#!/usr/bin/Rscript

#Module membership calculations
#---
#Script to calculate the module memberships of genes for each module obtained from a WGCNA network.
#Module membership is the correlation between gene expression values and module eigengenes.


#load packages
library(dplyr)
library(jsonlite)

#define global variables
DIR = "C:/Users/dstra/Nextcloud/Documents/Computers/Code/Github repos/bioinformatics analyses/Machine learning/Immune cell classification" #working directory

#set working directory
setwd(DIR)

#load data
network = fromJSON("Data/signed_module_network.json")
scaled_counts = read.csv("Data/TPM_scaled_counts.csv", row.names = 1)
metadata = read.csv("Data/metadata.csv")


#Define network results
#---
#define the module eigengene data
ME_data = network$MEs

#define the gene modules
modules = network$colors
names(modules) = network$genes

#define the module names
module_names = unique(modules)


#Correlate gene expression values with module eigengenes (calculate module memberships)
#---
#define an empty list
MM_results = list()
for (module_name in module_names){
	#define the gene list for the module
	module_genes = names(modules[modules %in% module_name])

	#define the VST gene expression values
	scaled_values = filter(scaled_counts, rownames(scaled_counts) %in% module_genes)

	#define the module MEs
	MEs = ME_data[, paste0("ME", module_name)]

	#correlate MEs with VST expression values
	MMs = apply(scaled_values, 1, function(row) cor(row, MEs, use = 'complete.obs', method = 'pearson'))

	#store the results
	MM_df = data.frame(ensembl = names(MMs), MM = MMs)
	rownames(MM_df) = NULL
	MM_results[[module_name]] = MM_df
}

#save results as JSON
write(toJSON(MM_results), "Data/module_memberships.json")
