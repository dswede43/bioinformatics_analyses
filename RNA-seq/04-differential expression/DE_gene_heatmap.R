#!/usr/bin/Rscript

#Differentially expressed gene LFC heatmap
#---
#This script creates a heatmap to represent the log2 Fold Changes (LFC) from differentially expressed (DE)
#genes across multiple conditions or comparisons.


#Table of contents:
#---
#1. Ensure that the lists of gene groups are differentially expressed
#2. Obtain the gene group LFC values
#3. Convert the LFC results into matrices
#4. Plot heatmaps

#load packages
library(dplyr)
library(tibble)
library(openxlsx)
library(jsonlite)
library(pheatmap)

#define global variables
DIR = "/path/to/directory" #working directory
COMPARISONS = c("differential","expression","comparisons") #define the differential expression comparisons
LFC_CUTOFF = 0.5 #log2 fold change cutoff
PADJ_CUTOFF = 0.05 #adjusted p-value cutoff
BREAKS = seq(-2, 2, by = 0.1)
COLORS = colorRampPalette(c("#FFDD00","#FFFFFF","#00BBFF"))(length(BREAKS))
HEIGHT = 8 #heatmap cell height
WIDTH = 12 #heatmap cell width
FONT = 8 #heatmap font size

#set working directory
setwd(DIR)

#load data
DE_results = fromJSON("DESeq2_results.JSON") #list of DE gene result data frames

gene_names = read.csv("gene_identities.csv") #gene annotations

gene_group_names = getSheetNames("gene_lists.xlsx") #lists of gene groups
gene_groups = list()
for (gene_group_name in gene_group_names){
	gene_groups[[gene_group_name]] = read.xlsx("genes_lists.xlsx", sheet = gene_group_name, colNames = FALSE)
}


#1. Ensure that the lists of gene groups are differentially expressed
#---
#define an empty array
DE_genes = c()
for (comparison in COMPARISONS){
	#define the LFC results for the current comparison
	LFC_result = DE_results[[comparison]]

	#filter for differentially expressed genes
	LFC_result = filter(LFC_result, abs(log2FoldChange) > LFC_CUTOFF & padj < PADJ_CUTOFF)

	#store the results
	DE_genes = c(DE_genes, LFC_result$ensembl)
}

#determine all unique genes
DE_genes = unique(DE_genes)

#convert gene identities
DE_genes = filter(gene_names, ensembl %in% DE_genes)$hgnc
DE_genes = DE_genes[DE_genes != ""]

#define an empty list
DE_gene_groups = list()
for (gene_group in names(gene_groups)){
	#define the genes for the current group
	genes = gene_groups[[gene_group]]$X1

	#retain genes that are differentially expressed
	DE_gene_groups[[gene_group]] = intersect(genes, DE_genes)
}


#2. Obtain the gene group LFC values
#---
#define an empty list
LFC_lists = list()
for (gene_group in names(DE_gene_groups)){
	#define the current gene group
	genes = DE_gene_groups[[gene_group]]

	#convert gene identities
	genes = filter(gene_names, hgnc %in% genes)$ensembl

	for (comparison in COMPARISONS){
		#define the LFC results for the current comparison
		LFC_result = DE_results[[comparison]]

		#filter for gene groups
		LFC_result = LFC_result %>%
			filter(ensembl %in% genes) %>%
			select(ensembl, log2FoldChange)

		#store the results
		LFC_lists[[gene_group]][[comparison]] = LFC_result
	}
}


#3. Convert the LFC results into matrices
#---
#define an empty list
LFC_mats = list()
for (gene_group in names(LFC_lists)){
	#define the current gene group
	genes = DE_gene_groups[[gene_group]]

	#define the current LFC values
	LFC_list = LFC_lists[[gene_group]]

	#convert to matrix
	LFC_mat = do.call(cbind, LFC_list)
	LFC_mat = LFC_mat[, grep("log2FoldChange", colnames(LFC_mat))]
	rownames(LFC_mat) = genes
	colnames(LFC_mat) = COMPARISONS

	#store the results
	LFC_mats[[gene_group]] = LFC_mat
}


#4. Plot heatmaps
#---
for (gene_group in names(LFC_mats)){
	#define the current LFC matrix
	LFC_mat = LFC_mats[[gene_group]]

	#plot the heatmap
	pheatmap(LFC_mat,
		 border_color = "NA",
		 cluster_cols = FALSE,
		 cluster_rows = FALSE,
		 angle_col = 45,
		 color = COLORS,
		 breaks = BREAKS,
		 cellheight = HEIGHT,
		 cellwidth = WIDTH,
		 fontsize = FONT,
		 filename = paste0("heatmap_", gene_group, ".pdf"))
}
