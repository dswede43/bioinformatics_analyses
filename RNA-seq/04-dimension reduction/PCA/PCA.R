#!/usr/bin/Rscript

#Principle Component Analysis (PCA)
#---
#Script to run a PCA on RNA-sequencing samples using the variance stabilized read counts.

#define global variables
dir = "path/to/directory"
norm_count_file = 'your_norm_read_counts.csv' #normalized read count matrix file name
metadata_file = 'your_metadata.csv' #metadata file name
vst_count_file = 'your_vst_read_counts.csv' #VST read count matrix file name
biotypes_file = 'your_gene_biotypes.csv' #gene biotypes
gene_names_file = 'your_gene_names.csv' #gene names

gene_biotype = 'protein_coding' #gene biotypes of interest
threshold = 5 #normalized read count exclusion threshold
variables = c('your','categorical','variables') #categorical variables

#set working directory
setwd(dir)

#load packages
library(dplyr)
library(ggplot2)
library(ggpubr)
source("https://raw.githubusercontent.com/dswede43/bioinformatics_analyses/main/Functions/load_read_count_data.R")

#load data
data = load_read_count_data(norm_count_file, metadata_file, variables)
norm_counts = data[['read_counts']] #normalized read count matrix
metadata = data[['metadata']] #sample metadata
vst_df = read.csv(vst_count_file, header = TRUE, sep = ",", row.names = 1) #variance stabilized transformed (VST) normalized read counts
biotypes = read.csv(biotypes_file, header = TRUE, sep = ",") #gene RNA biotypes
gene_names = read.csv(gene_names_file, header = TRUE, sep = ",")


#Remove genes based on low normalized read counts and biotypes
#---
#remove genes with average normalized read counts below the threshold
genes = rownames(norm_counts[rowMeans(norm_counts[]) > threshold,])

#remove uninteresting gene biotypes
genes = filter(biotypes, biotype == gene_biotype & ensembl %in% genes)[,'ensembl']

#subset the VST matrix for the genes of interest
vst_df = filter(vst_df, rownames(vst_df) %in% genes)


#Perform the PCA
#---
#calculate PC scores for each sample without scaling (VST transformation has already scaled the normalized read counts)
pca = prcomp(t(vst_df), scale = FALSE)

#obtain the variability explained by each PC
var = pca$sdev^2
per = round(var/sum(var)*100, 1)

#create a data frame of PC scores
pca_df = data.frame(metadata, PC1 = pca$x[ ,1], PC2 = pca$x[ ,2]) #create a PC data frame

#obtain PC loading scores
loading_scores = data.frame(pca$rotation[,1:2])


#Find the most influential genes
#---
#define a function to calculate the hypotenuse
calculate_hypotenuse = function(row){
	hypotenuse = sqrt(row[1]^2 + row[2]^2)
	return(hypotenuse)
}

#obtain the top most influential genes explaining variance
top_genes = names(sort(apply(loading_scores, 1, calculate_hypotenuse), decreasing = TRUE))[1:10]

#subset the loading scores to the most influential genes
top_loading_scores = filter(loading_scores, rownames(loading_scores) %in% top_genes)
top_loading_scores = data.frame(ensembl = rownames(top_loading_scores), top_loading_scores)
rownames(top_loading_scores) = NULL

#add gene names
top_loading_scores = merge(top_loading_scores, gene_names[,c('ensembl','hgnc')], by = 'ensembl')


#Visualize PCA bi-plots
#---
plot_list = list()

#colour the points of each PCA plot by the variables of interest
for(variable in variables){
    plot = ggplot() +
	
	#create scatterplot
    geom_point(data = pca_df, aes_string(x = 'PC1', y = 'PC2', col = variable),
		size = 2) +
	labs(
		x = paste('PC1 - ', per[1], '%', sep = ''),
		y = paste('PC2 - ', per[2], '%', sep = '')) +

	#add loading score vectors
	geom_segment(data = top_loading_scores,
		aes(x = 0, y = 0, xend = PC1 * 1000, yend = PC2 * 1000),
		arrow = arrow(length = unit(0.2, "cm")), size = 0.1, col = 'black') +
	
	#label each loading score
	geom_text(data = top_loading_scores,
		aes(x = PC1 * 1000, y = PC2 * 1000, label = hgnc),
		size = 3) +
    theme_classic()

	#store plot
    plot_list[[variable]] = plot
}

#visualize all variables in one plot
pca_plots = ggarrange(plotlist = plot_list)
pca_plots

#save plot as pdf
pdf("PCA_plots.pdf", height = 8, width = 8)
pca_plots
dev.off()
