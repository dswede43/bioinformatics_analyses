#!/usr/bin/Rscript

#Module hub genes
#---
#Script to identify the intra-modular hub genes based on module membership (MM) values obtained from WGCNA.
#Genes are considered a "hub gene" if their MM absolute value is greater than a specific percentile threshold.

#load packages
library(dplyr)
library(jsonlite)

#define global variables
DIR = "/path/to/directory" #working directory
MM_PERCENTILE = 0.9 #module membership (MM) percentile threshold

#set working directory
setwd(DIR)

#load data
MM_results = fromJSON("module_memberships.json")


#Identify intra-modular hub genes
#---
#define the current study phase
MM_values = MM_results[[phase]]

#define empty list
hub_genes = list()
for (module_name in names(MM_values)){
	#define the current module
	module_genes = MM_values[[module_name]]

	#identify module hub genes
	module_genes = module_genes %>%
		filter(abs(MM) > quantile(abs(MM), probs = MM_PERCENTILE, na.rm = TRUE)) %>%
		arrange(desc(MM))

	#store the results
	hub_genes[[module_name]] = module_genes
}

#save results as xlsx
write.csv(hub_genes, "module_hub_genes.csv")
