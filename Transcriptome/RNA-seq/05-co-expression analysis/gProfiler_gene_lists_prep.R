#!/usr/bin/Rscript

#Multiquery input preparation for gProfiler
#---
#Script to prepare the list of genes modules obtained from WGCNA for a multiquery input to
#g:Profiler g:GOSt. Gene ensembl id's are ordered by module membership to complete the
#functional enrichment of gene modules.

#load packages
library(dplyr)
library(jsonlite)

#define global variables
DIR = "/path/to/directory" #working directory
MM_CUTOFF = 0.1

#set working directory
setwd(DIR)

#load data
MM_results = fromJSON("module_memberships.json")


#Create list of genes for input with g:Profiler g:GOSt
#---
#define the module names
module_names = names(MM_results)

#define empty list
gene_lists = c()
for (module_name in module_names){
	#subset to the gene module
	module_genes = MM_results[[module_name]]

	#remove MMs less than the cutoff and order genes in decreasing order of MM
	module_genes = module_genes %>%
		filter(abs(MM) > MM_CUTOFF) %>%
		arrange(desc(MM))

	#append the list of ensembl id's to growing list
	gene_lists = c(gene_lists, c(paste0(">", module_name), module_genes$ensembl))
}

#store the results as TXT
write(gene_lists, "gProfiler_multiquery_input.txt")
