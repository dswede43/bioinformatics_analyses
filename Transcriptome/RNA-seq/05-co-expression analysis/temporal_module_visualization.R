#!/usr/bin/Rscript

#Module z-score plot
#---
#Script to visualize the expression changes of modules obtained from WGCNA for a temporal study design.
#The visualization is created by Z-score scaling the normalized read counts across genes at each timepoint.
#Z-scores are then plotted across time to visualize the expression changes for each gene module.

#Table of contents:
#---
#1. Filter to expressed protein-coding genes
#2. Scale the read count data
#3. Create the z-score line plot

#load packages
library(dplyr)
library(jsonlite)
library(ggplot2)

#define global variables
DIR = "/path/to/directory" #working directory
SAMPLES = 'all_samples' #sample names

#set working directory
setwd(DIR)

#load data
networks = fromJSON("module_networks.json")
counts = read.csv("norm_read_counts.csv", row.names = 1)
metadata = read.csv("metadata.csv")
gene_biotypes = read.csv("gene_biotypes.csv")


#1. Filter to expressed protein-coding genes
#---
#obtain list of protein-coding genes
protein_coding_genes = filter(gene_biotypes, biotype == 'protein_coding')$ensembl

#filter expression matrix for protein-coding genes
counts = filter(counts, rownames(counts) %in% protein_coding_genes)


#2. Scale the read count data
#---
#scale read counts
scale_counts = t(scale(t(counts)))
scale_counts = data.frame(ensembl = rownames(scale_counts), scale_counts)
rownames(scale_counts) = NULL

#convert data to long format
count_values = reshape2::melt(scale_counts, id = 'ensembl')
colnames(count_values) = c('ensembl','sample','count')
count_values = merge(count_values, metadata[, c(1,7)], by = 'sample')


#3. Create the z-score line plot
#---
#define the network
network = networks[[SAMPLES]]

#create data frame of modules
modules_df = data.frame(ensembl = network$genes, module = network$colors)

#summarize and average gene counts
modules_summary = count_values %>%
	group_by(ensembl, time) %>%
	summarise(mean_count = mean(count), nrep = n())

#join module information
modules_summary = full_join(modules_summary, modules_df, by = 'ensembl')

#define module colors
colors = unique(modules_summary$module)
modules_summary$module = factor(modules_summary$module, levels = colors)

#create z-score line plot
line_plot = ggplot(data = modules_summary, aes(x = time, y = mean_count)) +
	geom_line(aes(group = ensembl), linewidth = 0.2, color = 'grey20', alpha = 0.01) +
	geom_line(aes(group = module, color = module), linewidth = 1, stat = 'summary', fun = 'median') +
	labs(x = 'timepoints (days)', y = 'mean scaled expression (Z-score)') +
	scale_color_manual(values = colors) +
	theme_classic()

#save plot as pdf
pdf(paste0("zscore_plot-", SAMPLES, ".pdf"), height = 6, width = 10)
line_plot
dev.off()
