#!/usr/bin/Rscript

#Data distribution visualization
#---
#Script to visualize the RNA-seq data distributions.


#Table of contents:
#---
#1. Read count distribution by immune cell types
#2. Immune cell type differences with PCA


#load packages
library(dplyr)
library(reshape2)
library(ggplot2)
library(introdataviz)

#define global variables
DIR = "C:/Users/dstra/Nextcloud/Documents/Computers/Code/Github repos/bioinformatics analyses/Machine learning/Immune cell classification" #working directory
LOADING_SCORE_NUM = 5 #number of variables to display in PCA bi-plot
LOADING_SCORE_SIZE = 300 #loading score multiplier for plotting

#set working directory
setwd(DIR)

#load data
count_data = read.csv("Data/ML_train_data.csv") #training data
gene_names = read.csv("Data/gene_biotypes-names.csv")[-3] #gene names


#1. Read count distribution by immune cell types
#---
#format the data for plotting
plot_data = melt(count_data, id = 'cell_group')
colnames(plot_data) = c('cell_group','ensembl','read_count')
plot_data$cell_group = ifelse(plot_data$cell_group == 'T_cell', 'T-cell', 'non T-cell')

#plot the read count distributions by immune cell type
flat_violin = ggplot(data = plot_data, aes(x = "", y = log2(read_count), fill = cell_group)) +
	introdataviz::geom_flat_violin(trim = FALSE, alpha = 0.4,
		position = position_nudge(x = 0.5)) +
	geom_point(aes(color = cell_group), size = 0.5, alpha = 0.01, show.legend = FALSE,
		position = position_jitter(height = 0)) +
	geom_boxplot(width = 0.1, alpha = 0.5, show.legend = FALSE, outlier.shape = NA,
		position = position_nudge(x = -0.5)) +
	stat_summary(fun.data = mean_se, mapping = aes(color = cell_group), show.legend = FALSE,
		position = position_nudge(0.7)) +
	labs(y = 'log2 gene expression values (TPM normalized read counts)', x = NULL) +
	scale_fill_discrete(name = 'immune cell type') +
	coord_flip() +
	theme_classic() +
	theme(panel.grid.major.y = element_blank(),
		plot.margin = unit(c(1,1,2,2), 'cm'))

#save plot as PDF
pdf("Visualizations/gene_expression_distribution.pdf", height = 8, width = 8)
flat_violin
dev.off()


#2. Immune cell type differences with PCA
#---
#calculate PC scores
pca = prcomp(count_data[-1], scale = TRUE)

#obtain the variability explained by each PC
var = pca$sdev^2
per = round(var/sum(var)*100, 1)

#create a data frame of PC scores
pca_df = data.frame(cell_group = count_data$cell_group, PC1 = pca$x[ ,1], PC2 = pca$x[ ,2]) #create a PC data frame
pca_df$cell_group = ifelse(pca_df$cell_group == 'T_cell', 'T-cell', 'non T-cell')

#obtain PC loading scores
loading_scores = data.frame(pca$rotation[,1:2])

#define a function to calculate the hypotenuse
calculate_hypotenuse = function(row){
	hypotenuse = sqrt(row[1]^2 + row[2]^2)
	return(hypotenuse)
}

#obtain the top most influential genes explaining variance
top_genes = names(sort(apply(loading_scores, 1, calculate_hypotenuse), decreasing = TRUE))[1:LOADING_SCORE_NUM]

#subset the loading scores to the most influential genes
top_loading_scores = filter(loading_scores, rownames(loading_scores) %in% top_genes)
top_loading_scores = data.frame(ensembl = rownames(top_loading_scores), top_loading_scores)
rownames(top_loading_scores) = NULL

#add gene names
top_loading_scores = merge(top_loading_scores, gene_names, by = 'ensembl')

#visualize the PCA bi-plot
pca_plot = ggplot() +
	geom_point(data = pca_df, aes(x = PC1, y = PC2, color = cell_group), size = 2, alpha = 0.4) +
	geom_segment(data = top_loading_scores,
		aes(x = 0, y = 0, xend = PC1 * LOADING_SCORE_SIZE, yend = PC2 * LOADING_SCORE_SIZE),
		arrow = arrow(length = unit(0.1, "cm")), linewidth = 0.1, col = 'black', alpha = 0.6) +
	geom_text(data = top_loading_scores,
		aes(x = PC1 * LOADING_SCORE_SIZE, y = PC2 * LOADING_SCORE_SIZE, label = hgnc),
		size = 2, alpha = 0.8) +
	labs(x = paste0('PC1 - ', per[1], '%'), y = paste0('PC2 - ', per[2], '%')) +
	scale_color_discrete(name = 'immune cell type') +
	theme_classic()

#save plot as PDF
pdf("Visualizations/PCA_gene_expression.pdf", height = 5, width = 7)
pca_plot
dev.off()
