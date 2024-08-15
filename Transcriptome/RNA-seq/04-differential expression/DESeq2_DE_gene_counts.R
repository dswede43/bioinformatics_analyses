#!/usr/bin/Rscript

#Counting differentially expressed (DE) genes for pairwise comparisons from DESeq2 results
#---
#Script to query and count the number of differential expressed genes resulting from pairwise comparisons
#obtained with DESeq2.

#load packages
library(dplyr)
library(jsonlite)
library(openxlsx)

#define global variables
DIR = "/path/to/directory" #working directory
GENE_BIOTYPES = c("protein_coding","IG_C_gene","TR_C_gene","IG_V_gene","TR_V_gene") #gene biotypes to keep
COMPARISONS = c("pairwise","comparison","names")
LFC_CUTOFF = 0.5 #log2 fold change cutoff
PADJ_CUTOFF = 0.05 #adjusted p-value cutoff

#set working directory
setwd(DIR)

#load data
DE_results = fromJSON("DESeq2_table_results.JSON")
gene_biotypes = read.csv("RNA_biotypes.csv")
gene_names = read.csv("Gene_identities.csv")


#Summarize the number of differential expressed genes for each pairwise comparison
#---
#define an empty data frame
DE_counts = data.frame()
for (comparison in COMPARISONS){
	#define current pairwise comparison
	DE_result = DE_results[[comparison]]

	#filter for DE genes
	DE_result = filter(DE_result, abs(log2FoldChange) > LFC_CUTOFF & padj < PADJ_CUTOFF)

	#merge HGNC symbols
	DE_result = merge(DE_result, gene_names, by = 'ensembl')

	#merge gene biotypes
	DE_result = merge(DE_result, gene_biotypes, by = 'ensembl')

	#store the results
	DE_result = dplyr::select(DE_result, ensembl, hgnc, biotype, log2FoldChange, padj)

	#obtain the total number of DE genes
	n_genes = nrow(DE_result)

	#obtain the total number of up-regulated genes
	n_genes_up = nrow(filter(DE_result, log2FoldChange > 0))

	#obtain the total number of down-regulated genes
	n_genes_down = nrow(filter(DE_result, log2FoldChange < 0))

	#obtain the number of DE genes of each RNA biotype
	n_biotypes = c()
	for (gene_biotype in GENE_BIOTYPES){
		#obtain the number of DE genes for the current RNA biotype
		tmp = nrow(filter(DE_result, biotype %in% gene_biotype))

		#store the results
		n_biotypes = c(n_biotypes, tmp)
	}

	#append the RNA biotype names to the array
	names(n_biotypes) = GENE_BIOTYPES

	#obtain the total number of DE genes from all biotypes
	total_biotypes = sum(n_biotypes)

	#obtain the total number of up-regulated genes from all biotypes
	biotypes_up = nrow(filter(DE_result, biotype %in% GENE_BIOTYPES & log2FoldChange > 0))

	#obtain the total number of down-regulated genes from all biotypes
	biotypes_down = nrow(filter(DE_result, biotype %in% GENE_BIOTYPES & log2FoldChange < 0))

	#store the results
	DE_counts = rbind(DE_counts, data.frame(comparison,
											n_genes,
											n_genes_up,
											n_genes_down,
											t(n_biotypes),
											total_biotypes,
											biotypes_up,
											biotypes_down))
}

#save results as CSV
write.csv(DE_counts, "DE_gene_counts.csv", row.names = FALSE)
