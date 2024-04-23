#!/usr/bin/Rscript

#LFC GSEA dotplot
#---
#This script creates a dotplot to visualize the results from a
#Gene Set Enrichment Analysis (GSEA) using log2 fold changes (LFC).

#define global variables
dir = "/path/to/directory" #directory to save data
GSEA_results_file = 'your_GSEA_results.xlsx' #file name of GSEA results

comparisons = c('your','differential','expression','comparisons') #define differential expression (DE) comparisons (***must be the same as column header names***)
pos_NES_cutoff = 2 #cutoff for positve NES
neg_NES_cutoff = -2 #cutoff for negative NES
dotplot_cutoff = 5 #cutoff for the number of terms to display in the dotplot

#set working directory
setwd(dir)

#load packages
library(openxlsx)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggvenn)

#load data
gsea_results = list()
#GSEA result tables
for(comparison in comparisons){

	#load the GSEA results for a given DE comparison
	gsea_results[[comparison]] = read.xlsx(GSEA_results_file, sheet = comparison)
}


#Dot plot
#---
dotplot_list = list()

for(comparison in comparisons){

	#subset the GSEA results to a given DE comparison
	gsea_df = gsea_results[[comparison]]
	
	#if there exists GO terms
	if(nrow(gsea_df) != 0){
	
		#obtain the GO terms with the most positive and negative enrichment scores
		gsea_df = gsea_df %>%
			arrange(desc(NES)) %>%
			slice(c(1:dotplot_cutoff, (n() - (dotplot_cutoff - 1)):n())) %>%
			mutate(comparison = comparison)
		
		#store the results
		dotplot_list[[comparison]] = gsea_df
	}
}
#create a single data frame of results
dotplot_df = do.call(rbind, dotplot_list)
dotplot_df = mutate(dotplot_df, label = paste0(Description, " (", setSize, ")")) #add gene setSizes to the terms
rownames(dotplot_df) = NULL

#create the dotplot
dot_plot = ggplot(dotplot_df, aes(x = factor(comparison, levels = comparisons), y = factor(label, levels = unique(label)))) +
	geom_point(aes(color = NES, size = log10(p.adjust))) +
	labs(x = NULL, y = NULL) +
	scale_y_discrete(limits = rev) + #define y-axis order
	scale_color_gradientn(colours = c("red","white","green"), #define dot color gradient
		guide = guide_colorbar(order = 1)) +
	labs(color = "NES", size = "log10(p adjusted)") + #add legend titles
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) #slant the x-axis tick labels

#save plot as pdf
pdf("GSEA_dotplot.pdf", height = 8, width = 8)
dot_plot
dev.off()
