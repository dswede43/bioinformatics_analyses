#!/usr/bin/Rscript

#WGCNA phenotypic correlations
#---
#Script to identify modules with significant phenotypes through the correlation of clinical
#variables with module eigengenes.


#Table of contents:
#---
#1. Clean and format the clinical variable dataset
#2. Check for normality assumption among ME and clinical variables
#3. Scale and correlate module eigengene values to clinical variables
#4. Calculate correlation p-values (student asymptotic p-value)
#5. Create heatmaps to visualize correlations

#load packages
library(dplyr)
library(reshape2)
library(jsonlite)
library(ggplot2)
library(WGCNA)

#define global variables
DIR = "/path/to/directory" #working directory

CLIN_VARS = c("various", #clinical variables
			  "clinical",
			  "variables")

#set working directory
setwd(DIR)

#load data
clinical_data = read.csv("clinical_data.csv", strip.white = TRUE) #clinical variable data
network = fromJSON("module_network.json")


#1. Clean and format the clinical variable dataset
#---
#define columns to remove
clinical_data = filter(clinical_data, parameter %in% CLIN_VARS)
clinical_data = type.convert(clinical_data, as.is = TRUE)
clinical_data = reshape2::melt(clinical_data)


#2. Check for normalilty assumption among ME and clinical variables
#---
#define the ME data
ME_data = network$MEs

#convert data into long format
ME_data_long = reshape2::melt(ME_data)

#create the histogram
ME_histos = ggplot(data = ME_data_long, aes(x = value)) +
	geom_histogram(bins = 30) +
	labs(x = "module eigengenes (ME)", title = phase) +
	theme_classic()

#save plots as pdf
pdf("ME_histogram.pdf", height = 6, width = 6)
 ME_histos
dev.off()

#clinical variables
clin_histos = list()
for (clin_var in CLIN_VARS){
	#subset to current clinical variable
	clin_var_df = filter(clinical_data, parameter == clin_var)

	#create the histogram
	clin_histos[[clin_var]] = ggplot(data = clin_var_df, aes(x = measure)) +
		geom_histogram(bins = 30) +
		labs(x = clin_var) +
		theme_classic()
}

#save plots as pdf
pdf("Clinical_vars_histograms.pdf", height = 8, width = 8)
ggarrange(plotlist = clin_histos, ncol = 3, nrow = 5)
dev.off()


#3. Scale and correlate module eigengene values to clinical variables
#---
#define empty list
corr_results = list()

#scale the module eigengene data
ME_data = scale(ME_data)

for (clin_var in CLIN_VARS){
	#define and scale the clinical variable
	clin_var_df = filter(clinical_data, parameter == clin_var)
	clin_var_df = data.frame(clin_var_df, scaled_measure = scale(clin_var_df$measure))

	if (all(clin_var_df$measure == 0) | all(is.na(clin_var_df$scaled_measure))){
		#define correlations as NA
		corr_results[[clin_var]] = rep(NA, dim(ME_data)[2])
	}

	else {
		#correlate MEs with clinical variable
		corr_results[[clin_var]] = apply(ME_data, 2, function(col) cor(col, clin_var_df$scaled_measure, use = 'complete.obs', method = 'pearson'))
	}
}


#4. Calculate correlation p-values (student asymptotic p-value)
#---
#define empty list
pval_results = list()

for(clin_var in CLIN_VARS){
	#calculate the correlation p-value
	pval_results[[clin_var]] = corPvalueStudent(cor = corr_results[[clin_var]], nSamples = dim(ME_data)[1])
}


#5. Create heatmaps to visualize correlations
#---
#define correlation matrix
corr_mat = do.call(cbind, corr_results)

#define p-value matrix
pval_mat = do.call(cbind, pval_results)
pval_mat[is.na(pval_mat)] = 0 #replace NA values with zero

#add asterisks to correlation values
asterisks_mat = round(corr_mat, 2)
asterisks_mat[pval_mat < 0.05] = paste0(asterisks_mat[pval_mat < 0.05], "*")
asterisks_mat[asterisks_mat == "NA*"] = "NA"

#prepare data for plotting
heatmap_mat = melt(corr_mat)
heatmap_mat = data.frame(heatmap_mat, melt(asterisks_mat)$value)
colnames(heatmap_mat) = c("module","clin_var","correlation","asterisks")

#plot the heatmap
heatmap = ggplot(data = heatmap_mat, aes(x = clin_var, y = module, fill = correlation)) +
	geom_tile() +
	geom_text(aes(label = asterisks), size = 2) +
	scale_fill_gradient2(low = "blue1", high = "red", mid = "white", 
		midpoint = 0, limits = c(-1, 1), 
		breaks = seq(-1, 1, by = 0.5), 
		name = "Correlation") +
	labs(x = NULL, y = NULL, title = phase) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

#save plots as pdf
pdf("ME_corr_heatmap.pdf", height = 8, width = 8)
heatmap
dev.off()
