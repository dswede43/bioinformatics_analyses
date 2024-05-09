#!/usr/bin/Rscript

#Weighted Gene Co-expression Network Analysis (WGCNA)
#---
#Script to run a WGCNA for a bulk RNA-seq experiment. This script assumes the read count data has
#been normalized and transformed using the variance stabilizing transformation (VST). It also assumes
#genes with low read counts have already been filtered out.

#load packages
library(dplyr)
library(WGCNA)
library(ggplot2)
library(ggpubr)
library(jsonlite)

#define global variables
DIR = "/path/to/directory"
SIGN = 'signed'

#set working directory
setwd(DIR)

#load data
vst_counts = read.csv("VST_read_counts.csv", row.names = 1)


#Optimize the soft threshold power for scale-free netwotk topology
#---
#define powers to test
powers = c(1:20)

#optimize the soft threshold powers
exp_mat = t(vst_counts)

#calculate network topologies for each power
soft_thresholds = pickSoftThreshold(exp_mat, #expression matrix (rows = samples, columns = genes)
									powerVector = powers, #vector of potential soft power thresholds
									networkType = SIGN, #network type (signed or unsigned)
									verbose = 5)
soft_thresholds_df = soft_thresholds$fitIndices

#visualize the most optimial soft threshold power

#R-squared model fit
r_squared = ggplot(data = soft_thresholds_df[[phase]], aes(x = Power, y = SFT.R.sq, label = Power)) +
	geom_point(size = 2) +
	geom_text(nudge_y = 0.05) +
	geom_hline(yintercept = 0.8, color = 'red', linetype = 'dashed') +
	labs(title = phase, x = "powers", y = "scale-free fit (R-squared)") +
	theme_classic()

#mean connectivity
mean_connectivity = ggplot(data = soft_thresholds_df[[phase]], aes(x = Power, y = mean.k., label = Power)) +
	geom_point(size = 2) +
	labs(title = phase, x = "powers", y = "mean connectivity") +
	theme_classic()

#visualize the plots
plots = ggarrange(r_squared, mean_connectivity, ncol = 2)

#save the plots as pdf
pdf(paste0(SIGN, "_pickSoftThreshold.pdf"), height = 6, width = 6)
plots
dev.off()


#Identify modules using blockwiseModules funtion
#define the soft power threshold
soft_power = 10

#deal with function conflicts (use the WGCNA correlation function)
tmp_cor = cor
cor = WGCNA::cor

#construct the network in a block-wise manner
bw_network = blockwiseModules(exp_mat, #expression matrix (rows = samples, columns = genes)
							  maxBlockSize = 15000, #max number of genes for each block
							  TOMType = SIGN, #network type (signed or unsigned)
							  power = soft_power, #soft power to create weighted adjency matrix
							  mergeCutHeight = 0.25, #threshold to merge similar modules
							  numericLabels = FALSE, #module labels (TRUE = numbers, FALSE = colors)
							  deepSplit = 2, #sensitivity to module splitting (0 least and 4 most)
							  minModuleSize = 30, #min number of genes for a module
							  randomSeed = 123, #seed for RNG (for reproducible results)
							  verbose = 3)

#print the number of genes assigned to each module
table(bw_network$colors)

#save results as json
bw_network$genes = names(bw_network$colors)
write(toJSON(bw_networks), paste0(SIGN, "_module_network.json"))

#reassign correlation function
cor = tmp_cor
