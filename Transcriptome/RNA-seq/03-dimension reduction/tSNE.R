#!/usr/bin/Rscript

#t-Distributed Stochastic Neighbour Embedding
#---
#Script to run a t-SNE analysis for RNA-sequencing samples.

#Rtsne documentation:
#https://www.rdocumentation.org/packages/Rtsne/versions/0.17/topics/Rtsne


#Table of contents:
#---
#1. Run the t-SNE analysis on samples
#2. Visualize the t-SNE results

#load packages
library(Rtsne)
library(ggplot2)

#define global variables
DIR = "/path/to/directory" #working directory

#set working directory
setwd(DIR)

#load data
vst_counts = read.csv("VST_read_counts.csv", row.names = 1)
metadata = read.csv("metadata.csv")


#1. Run the t-SNE analysis on samples
#---
#set the seed to ensure reproducible results (t-SNE is stochastic)
set.seed(123)

#use the t-SNE algorithm on samples
tsne = Rtsne(as.matrix(t(vst_counts)),
			 dims = 2, #number of dimensions to output
			 perplexity = 12) #perplexity parameter (number of nearest neighbours consider)

#create the t-SNE data frame
tsne_df = data.frame(meta, tSNE1 = tsne$Y[ ,1], tSNE2 = tsne$Y[ ,2])


#2. Visualize the t-SNE results
#---
#create scatter plot
tsne_plot = ggplot(data = tsne_df, aes(x = tSNE1, y = tSNE2)) +
	geom_point(size = 3) +
	geom_text(vjust = 1.6, show.legend = F) +
	labs(color = "category",
		x = "t-SNE1",
		y = "t-SNE2") +
	theme_classic()

#save plot as PDF
pdf("tnse_plot.pdf", height = 6, width = 6)
tsne_plot
dev.off()
