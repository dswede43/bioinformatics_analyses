#!/usr/bin/Rscript

#Module summary
#---
#Script to summarize the modules obtained from WGCNA. The summary includes the module genes and their
#module membership (MM) values, hub genes, enrichment results, and module eigengene (ME) values.


#Table of contents:
#---
#1. Hub genes
#2. Enrichment results
#3. Module genes and module memberships (MM)
#4. Module eigengene (ME) values

#load packages
library(dplyr)
library(jsonlite)
library(openxlsx)

#define global variables
DIR = "/path/to/directory" #working directory
MM_CUTOFF = 0.1
MM_PERCENTILE = 0.9

#set working directory
setwd(DIR)

#load data
network = fromJSON("module_network.json")
module_memberships = fromJSON("module_memberships.json")
module_enrichment = fromJSON("module_enrichment_results.json")
gene_biotypes = read.csv("RNA_biotypes.csv")
gene_names = read.csv("gene_identities.csv")


#Function to summarize WGCNA results
#---
summarize_modules = function(network, module_memberships, module_enrichment, gene_biotypes, gene_names){
	#define an empty excel workbook
	wb = createWorkbook()

	#define an empty list
	module_summary = list()
	for (module_name in unique(network$colors)){

		#1. Hub genes
		#---
		#obtain module gene identities
		module_genes = module_memberships[[module_name]]

		#obtain the module hub genes
		hub_genes = module_genes %>%
			filter(abs(MM) > quantile(abs(MM), probs = MM_PERCENTILE, na.rm = TRUE)) %>%
			arrange(desc(MM))

		#store results
		module_summary[['hub_genes']] = hub_genes


		#2. Enrichment results
		#---
		#obtain the module enrichment results
		module_enr = module_enrichment[[module_name]]

		#store results
		module_summary[['module_enrichment']] = module_enr


		#3. Module genes and module memberships (MM)
		#---
		#remove genes below the MM cutoff
		module_genes = filter(module_genes, abs(MM) > MM_CUTOFF)

		#append gene biotype and HGNC symbols
		module_genes = merge(module_genes, gene_names, by = 'ensembl')
		module_genes = merge(module_genes, gene_biotypes, by = 'ensembl')

		#arrange genes in descending order of MM
		module_genes = arrange(module_genes, desc(MM))

		#store results
		module_summary[['module_genes']] = module_genes


		#4. Module eigengene (ME) values
		#---
		#define the module eigengene data
		ME_data = network$MEs

		#obtain module eigengene values
		MEs = dplyr::select(ME_data, contains(paste0("ME", module_name)))
		MEs = data.frame(sample = rownames(MEs), ME = MEs)
		rownames(MEs) = NULL

		#store results
		module_summary[['module_eigengenes']] = MEs

		#add a new sheet the excel workbook
		addWorksheet(wb, module_name)

		#define the starting rows for each data frame in the current sheet
		start_row = 1
		for (i in 2:(length(names(module_summary)) + 1)){
			start_row = c(start_row, ifelse(is.null(nrow(module_summary[[i - 1]])), 0, nrow(module_summary[[i - 1]])) + start_row[i - 1] + 3)
		}

		#define the labels for each data frame
		labels = names(module_summary)

		#write the data frames to the excel worksheet
		for (i in seq_along(module_summary)) {
			writeData(wb, sheet = module_name, x = labels[i], startRow = start_row[i], startCol = 1)
			writeData(wb, sheet = module_name, x = module_summary[[i]], startRow = start_row[i] + 1, startCol = 1)
		}
	}

	#save results as xlsx
	saveWorkbook(wb, file = "WGCNA_results_summary.xlsx"), overwrite = TRUE)
}


#Apply WGCNA summary function
#---
summarize_modules(network, module_memberships, module_enrichment, gene_biotypes, gene_names)
