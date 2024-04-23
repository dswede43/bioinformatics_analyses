#!/usr/bin/Rscript

#Differentially expressed gene LFC heatmap
#---
#This script creates a heatmap to represent the log2 Fold Changes (LFC) from differentially expressed (DE)
#genes across multiple conditions or comparisons.

#load packages
library(dplyr)
library(tibble)
library(jsonlite)
library(pheatmap)

#define global variables
dir = "path/to/directory" #working directory
gene_list_file = "gene_list_file.txt" #file name of gene list
gene_names_file = "gene_names_file.csv" #file name of gene HGNC symbols
DE_results_file = "DE_results_file.JSON" #file name of DE results

conditions = c("experimental","conditions") #study phase experimental conditions

#set working directory
setwd(dir)

#load data
gene_list = read.delim(gene_list_file, header = FALSE)$V1 #list of HGNC DE gene symbols
gene_names = read.csv(gene_names_file) #data frame of gene name conversions
DE_results = fromJSON(DE_results_file) #list of DE gene result data frames


#Create LFC matrix
#---
#define the differential expression comparisons
comparisons = names(DE_results)

#convert gene list from HGNC symbols to ensembl id's
gene_ensembl = filter(gene_names, hgnc %in% gene_list)$ensembl

LFC_list = list()
for (comparison in comparisons){
	#define the data frame for the current DE comparison
	LFC_result = DE_results[[comparison]]

	#filter for genes in the gene list
	LFC_result = filter(LFC_result, ensembl %in% gene_ensembl)

	#extract the LFC values of the genes
	LFC_list[[comparison]] = LFC_result[,c("ensembl","log2FoldChange")]
}

#convert into LFC results into matrix format
LFC_matrix = do.call(cbind, LFC_list)
rownames(LFC_matrix) = LFC_matrix[,1]
LFC_matrix = LFC_matrix[, grep("log2FoldChange", colnames(LFC_matrix))]
colnames(LFC_matrix) = comparisons

#convert gene names back to HGNC symbols
LFC_matrix = rownames_to_column(LFC_matrix, var = "ensembl")
LFC_matrix = left_join(LFC_matrix, gene_names[,c(1,3)], by = "ensembl")
LFC_matrix = column_to_rownames(LFC_matrix, var = "hgnc")[,-1]


#Create heatmap using pheatmap
#---
#define the heatmap aesthetics
breaks <- seq(min(LFC_matrix), max(LFC_matrix), by = 0.1) #set color scale breaks for heatmap
colors <- colorRampPalette(c("#FFDD00","#FFFFFF","#00BBFF"))(length(breaks)) #set color scale for heatmap
height <- 2 #set cellheight for heatmap
width <- 12 #set cellwidth for heatmap
font <- 8 #set fontsize for heatmap

#heatmap column annotations
annot_col <- data.frame(condition = rep(conditions, c(6, 7)))
rownames(annot_col) <- colnames(LFC_matrix)

#plot the heatmap
pheatmap(LFC_matrix,
		 color = colors,
		 border_color = "NA",
		 cluster_cols = FALSE,
		 cluster_rows = TRUE,
		 annotation_col = annot_col,
		 angle_col = 45,
		 breaks = breaks,
		 cellheight = height,
		 cellwidth = width,
		 fontsize = font)