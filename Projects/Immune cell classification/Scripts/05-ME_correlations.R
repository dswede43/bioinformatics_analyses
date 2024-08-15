#!/usr/bin/Rscript

#Module eigengene (ME) correlations with traits
#---
#Script to correlate WGCNA modules with phenotypic traits.


#Table of contents:
#---
#1. Prepate data
#2. Correlate modules with traits
#3. Visualize correlations
#4. Intramodular analysis: identify hub genes


#load packages
library(dplyr)
library(reshape2)
library(fastDummies)
library(WGCNA)
library(jsonlite)
library(ggplot2)

#define global variables
DIR = "C:/Users/dstra/Nextcloud/Documents/Computers/Code/Github repos/bioinformatics analyses/Machine learning/Immune cell classification" #working directory
PADJ_CUTOFF = 0.001 #adjusted p-value cutoff
MM_PERCENTILE = 0.9 #module membership (MM) percentile threshold

#set working directory
setwd(DIR)

#load data
network = fromJSON("Data/signed_module_network.json") #WGCNA modules
MM_data = fromJSON("Data/module_memberships.json") #module memberships
metadata = read.csv("Data/Metadata.csv") #study metadata
metadata = metadata[,c('sample','cell_group')]


#1. Prepate data
#---
#one-hot encode the traits
traits = dummy_cols(metadata,
					select_columns = "cell_group",
					remove_first_dummy = FALSE,
					remove_selected_columns = TRUE)

#format column and row names
colnames(traits) = gsub("cell_group_", "", colnames(traits))
rownames(traits) = traits$sample
traits = traits[,-1]

#module eigengene (ME) data
ME_data = network$MEs
colnames(ME_data) = gsub("ME", "", colnames(ME_data))


#2. Correlate modules with traits
#---
#define the number of sample 
n_samples = nrow(ME_data)

#correlate modules with traits
cor_mat = cor(ME_data, traits, use = 'p')

#calculate p-values for correlations
pval_mat = corPvalueStudent(cor_mat, n_samples)

#adjust for multiple comparisons
padj_mat = p.adjust(as.vector(as.matrix(pval_mat)), method = 'fdr')
padj_mat = matrix(padj_mat, nrow = nrow(pval_mat), ncol = ncol(pval_mat))
rownames(padj_mat) = rownames(pval_mat)
colnames(padj_mat) = colnames(pval_mat)


#3. Visualize correlations
#---
#round the correlations
annotated_mat = round(cor_mat, 2)

#add asterisks to correlation values
for (i in 1:nrow(annotated_mat)){
	for(j in 1:ncol(annotated_mat)){
		if (padj_mat[i, j] < 0.05){
			annotated_mat[i, j] = paste0(round(cor_mat[i, j], 2), "*")
		}
		if (padj_mat[i, j] < 0.01){
			annotated_mat[i, j] = paste0(round(cor_mat[i, j], 2), "**")
		}
		if (padj_mat[i, j] < 0.001){
			annotated_mat[i, j] = paste0(round(cor_mat[i, j], 2), "***")
		}	
	}
}

#prepare data for plotting
heatmap_mat = melt(cor_mat)
heatmap_mat = data.frame(heatmap_mat, melt(annotated_mat)$value)
colnames(heatmap_mat) = c("module","trait","correlation","asterisks")

#plot the heatmap
ME_trait_heatmap = ggplot(data = heatmap_mat, aes(x = trait, y = module, fill = correlation)) +
	geom_tile() +
	geom_text(aes(label = asterisks), size = 2) +
	scale_fill_gradient2(low = "blue1", high = "red", mid = "white", 
		midpoint = 0, limits = c(-1, 1), 
		breaks = seq(-1, 1, by = 0.5), 
		name = "Correlation") +
	labs(x = NULL, y = NULL) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

#save plot as PDF
pdf("Visualizations/ME_trait_cor_heatmap.pdf", height = 6, width = 8)
ME_trait_heatmap
dev.off()


#4. Intramodular analysis: identify hub genes
#---
padj_df = melt(padj_mat)
colnames(padj_df) = c("module","cell_group","padj")

sig_modules = as.character(filter(padj_df, cell_group == 'T_cell' & padj < PADJ_CUTOFF)$module)

#extract hub genes from each significant module
#define an empty list
hub_genes = list()
for (sig_module in sig_modules){
	#define the current module
	module_genes = MM_data[[sig_module]]

	#identify module hub genes
	module_genes = module_genes %>%
		filter(abs(MM) > quantile(abs(MM), probs = MM_PERCENTILE, na.rm = TRUE)) %>%
		arrange(desc(MM))

	#store the results
	hub_genes[[sig_module]] = module_genes
}

#save results as JSON
write(toJSON(hub_genes), "Data/module_hub_genes.json")
