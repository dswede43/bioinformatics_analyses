#!/usr/bin/Rscript

#Weighted Gene Co-expression Network Analysis (WGCNA)
#---
#Script to create gene modules using WGCNA.


#Table of contents:
#---
#1. Determine the most optimal soft thresholding power for scale-free network topology
#2. Identify modules using blockwiseModules() function
	#creates weighted adjency matrix using soft threshold power, convert to TOM matrix,
	#hierarchical clustering of dissimilarity matrix (1-TOM), cutreeDynamic() to determine modules

#load packages
library(dplyr)
library(WGCNA)
library(ggpubr)
library(ggplot2)
library(jsonlite)

#define global variables
DIR = "C:/Users/dstra/Nextcloud/Documents/Computers/Code/Github repos/bioinformatics analyses/Machine learning/Immune cell classification" #working directory
SIGN = 'signed' #network type for adjency matrix

#set working directory
setwd(DIR)

#load data
scaled_counts = read.csv("Data/TPM_scaled_counts.csv", row.names = 1)


#1. Optimize the soft threshold power for scale-free netwotk topology
#---
#define powers to test
powers = c(1:20)

#define and transpose the expression matrix
exp_mat = t(scaled_counts)

#calculate network topologies for each power
soft_thresholds = pickSoftThreshold(exp_mat, #expression matrix (rows = samples, columns = genes)
									powerVector = powers, #vector of potential soft power thresholds
									networkType = SIGN, #network type (signed or unsigned)
									verbose = 5)
soft_thresholds_df = soft_thresholds$fitIndices

#define the soft threshold power
soft_power = 12

#visualize the most optimial soft threshold power
#R-squared model fit
r_squared = ggplot(data = soft_thresholds_df, aes(x = Power, y = SFT.R.sq, label = Power)) +
	geom_point(size = 2) +
	geom_text(nudge_y = 0.05) +
	geom_hline(yintercept = 0.8, color = 'red', linetype = 'dashed') +
	geom_vline(xintercept = soft_power, color = 'blue', linetype = 'solid') +
	labs(x = "powers", y = "scale-free fit (R-squared)") +
	theme_classic()

#mean connectivity
mean_connectivity = ggplot(data = soft_thresholds_df, aes(x = Power, y = mean.k., label = Power)) +
	geom_point(size = 2) +
	geom_vline(xintercept = soft_power, color = 'blue', linetype = 'solid') +
	labs(x = "powers", y = "mean connectivity") +
	theme_classic()

#save the plots as pdf
pdf(paste0("Visualizations/", SIGN, "_pickSoftThreshold.pdf"), height = 4, width = 8)
ggarrange(r_squared, mean_connectivity, ncol = 2)
dev.off()


#2. Identify modules using blockwiseModules funtion
#---
#deal with function conflicts (use the WGCNA correlation function)
tmp_cor = cor
cor = WGCNA::cor

#construct the network in a block-wise manner
bw_network = blockwiseModules(exp_mat, #expression matrix (rows = samples, columns = genes)
							  maxBlockSize = 11000, #max number of genes for each block
							  TOMType = SIGN, #network type (signed or unsigned)
							  power = soft_power, #soft power to create weighted adjency matrix
							  mergeCutHeight = 0.25, #threshold to merge similar modules
							  numericLabels = FALSE, #module labels (TRUE = numbers, FALSE = colors)
							  deepSplit = 2, #sensitivity to module splitting (0 least and 4 most)
							  minModuleSize = 30, #min number of genes for a module
							  randomSeed = 123, #seed for RNG (for reproducible results)
							  #loadTOM = TRUE,
							  saveTOMs = TRUE, #save TOM adjency matrices
							  saveTOMFileBase = "Data/WGCNA_TPM_matrix", #TOM file name
							  verbose = 3)

#print the number of genes assigned to each module
table(bw_network$colors)

#visualize the modules with dendrogram and save as PDF
pdf("Visualizations/module_dendrogram.pdf", height = 6, width = 12)
plotDendroAndColors(dendro = bw_network$dendrograms[[1]],
					colors = bw_network$colors,
					dendroLabels = FALSE,
					addGuide = TRUE,
					autoColorHeight = FALSE,
					colorHeight = 0.1,
					hang = 0.03,
					guideHang = 0.05)
dev.off()

#store the results
bw_network$dendrograms = NULL
bw_network$genes = names(bw_network$colors)

#reassign correlation function
cor = tmp_cor

#save results as json
write(toJSON(bw_network), paste0("Data/", SIGN, "_module_network.json"))
