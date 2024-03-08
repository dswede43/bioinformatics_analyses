#!/usr/bin/Rscript

#Semantic similarity redundancy
#---
#This script reduces GO term redundancy using semantic similarity measures. It assumes an input of
#GO term identities with their normalized enrichment scores (NES) from GSEA. This method is inspired
#by the algorithm used in the REVIGO pipeline (https://doi.org/10.1371/journal.pone.0021800).

#define global variables
dir = "/path/to/directory"
GO_results_file = 'your_GO_results.csv' #name of GO results file

#set working directory
setwd(dir)

#load packages
library(GOSim)
library(dplyr)
library(ggplot2)
source('https://raw.githubusercontent.com/dswede43/bioinformatics_analyses/main/Functions/reduce_GO_list.R')

#load data
GO_results = read.csv(GO_results_file)


#Remove redundant GO terms
#---
#format data for semantic similarity reduction
GO_enr = GO_results$NES
names(GO_enr) = GO_results$ID

#define list of semantic similarity cutoffs
Cs = c(0.5,0.6,0.7,0.8,0.9,1)

GO_reduced_df = data.frame()
#for each cutoff
for (C in Cs){

	#count the number of GO terms after redundancy reduction
	GO_reduced_df = rbind(GO_reduced_df, c(C, length(reduce_GO_list(GO_enr, C = C))))
}
colnames(GO_reduced_df) = c('C', 'GO_count')


#Visualize results
#---
plot = ggplot(data = GO_reduced_df, aes(x = C, y = GO_count)) +
	geom_point(size = 3) +
	labs(x = "semantic similarity cutoff (C)", y = "number of GO terms") +
	theme_classic()
plot

#save plot as pdf
pdf("remove_redundant_GO_terms.pdf", height = 4, width = 4)
plot
dev.off()
