#!/usr/bin/Rscript

#WGCNA network module comparisons
#---
#Script to compare the gene modules obtained from WGCNA between two different networks.


#Table of contents:
#---
#1. Create the matrix of overlapping gene modules
#2. Calculate the statistical significance of the overlap (Fisher Exact Test)
#3. Visualize the module overlap matrix
#4. Summarize the gene identities

#load packages
library(dplyr)
library(jsonlite)
library(reshape2)
library(ggplot2)
library(openxlsx)

#define global variables
DIR = "/path/to/directory" #working directory
MM_CUTOFF = 0.1 #module membership cutoff
PADJ_CUTOFF = 0.05 #adjusted p-value cutoff
REF_LEN = 14807 #length of reference gene list

#set working directory
setwd(DIR)

#load data
networkA_MM = fromJSON("networkA-module_memberships.json")
networkB_MM = fromJSON("networkB-module_memberships.json")
gene_names = read.csv("gene_names.csv")


#1. Create the matrix of overlapping gene modules
#---
#define the module names from both studies
netB_modules = names(networkB_MM)
netA_modules = c(networkB_modules, setdiff(names(networkA_MM), netB_modules))

#define empty data frames
overlap_mat = data.frame()
pval_mat = data.frame()

#define empty arrays
netA_labels = c()
netB_labels = c()

#gene modules from network A
for (netA_module in netA_modules){
	#obtain the genes from the current module
	netA_genes = networkA_MM[[netA_module]]

	#remove genes below the module membership cutoff
	netA_genes = filter(netA_genes, abs(MM) > MM_CUTOFF)$ensembl

	#obtain the number of genes
	netA_len = length(netA_genes)

	#create the label for the current module
	netA_labels = c(netA_labels, paste0(netA_module, " (", netA_len, ")"))

	#define empty lists
	overlap_lens = c()
	pvals = c()

	#gene modules from network B
	for (netA_module in netB_modules){
		#obtain the genes from the current module
		netB_genes = networkB_MM[[netB_module]]

		#remove genes below the module membership cutoff
		netB_genes = filter(netB_genes, abs(MM) > MM_CUTOFF)$ensembl

		#obtain the number of genes
		netB_len = length(netB_genes)

		#create the label for the current module
		netB_labels = c(netB_labels, paste0(netB_module, " (", netB_len, ")"))

		#obtain the number of overlapping genes between both network modules
		overlap_len = length(intersect(netA_genes, netB_genes))

		#store the results
		overlap_lens = c(overlap_lens, overlap_len)


		#2. Calculate the statistical significance of the overlap (Fisher Exact Test)
		#---
		#create the contigency matrix
		mod = c(overlap_len, netB_len - overlap_len)
		tmp = netA_len - overlap_len
		non_mod = c(tmp, (REF_LEN - netB_len) - tmp)
		mat = matrix(c(mod, non_mod), 2,2, byrow = T)

		#test for statistical signficance using the Fisher Exact Test
		fisher_test = fisher.test(mat, alternative = "greater")
		pvals = c(pvals, fisher_test$p)
	}
	#define the module labels for study B
	netB_labels = unique(netB_labels)

	#store the results
	overlap_mat = rbind(overlap_mat, overlap_lens)
	pval_mat = rbind(pval_mat, pvals)
}

#adjust for multiple comparisons
padj_mat = p.adjust(as.vector(as.matrix(pval_mat)), method = 'fdr')
padj_mat = matrix(padj_mat, nrow = nrow(pval_mat), ncol = ncol(pval_mat))


#3. Visualize the module overlap matrix
#---
#add labels to overlap matrix
overlap_mat = data.frame(module = netA_labels, overlap_mat)
rownames(overlap_mat) = NULL
colnames(overlap_mat) = c('module', netB_labels)

#add labels to p-value matrix
padj_mat = data.frame(module = netA_labels, padj_mat)
rownames(padj_mat) = NULL
colnames(padj_mat) = c('module', netB_labels)

#transform data matrices into long format
overlap_long = melt(overlap_mat, id = 'module')
padj_long = melt(padj_mat, id = 'module')

#create the data frame for plotting
heatmap_df = cbind(overlap_long, padj_long[,3])
colnames(heatmap_df) = c('networkA','networkB','gene_count','pval')

#define the order of the heatmap cells
heatmap_df$networkA = factor(heatmap_df$networkA, levels = netA_labels)
heatmap_df$networkB = factor(heatmap_df$networkB, levels = rev(netB_labels))

#plot the heatmap
heatmap = ggplot(data = heatmap_df, aes(y = networkA, x = networkB, fill = -log(pval))) +
	geom_tile() +
	geom_text(aes(label = gene_count), color = 'black', size = 3) +
	scale_fill_gradient(low = 'white', high = 'red') +
	labs(y = 'networkA modules', x = 'networkB modules', title = 'networkA x networkB') +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

#save as pdf
pdf("Module_comparison_heatmap.pdf", height = 6, width = 8)
heatmap
dev.off()


#4. Summarize the gene identities
#---
#label the overlap and p-value matrices
rownames(overlap_mat) = overlap_mat$module
overlap_mat = overlap_mat[-1]
rownames(padj_mat) = padj_mat$module
padj_mat = padj_mat[-1]

#define an empty excel workbook
wb = createWorkbook()

#define an empty list
module_overlaps = list()

#gene modules from network A
for (i in 1:nrow(overlap_mat)){
	#define the current module
	netA_module = netA_modules[i]

	#obtain the genes from the current module
	netA_genes = networkA_MM[[netA_module]]

	#remove genes below the module membership cutoff
	netA_genes = filter(netA_genes, abs(MM) > MM_CUTOFF)$ensembl

	#define an empty list
	sig_overlaps = list()

	#gene modules from network B
	for (j in 1:ncol(overlap_mat)){
		#define the current module
		netB_module = netB_modules[j]

		#obtain the genes from the current module
		netB_genes = networkB_MM[[netB_module]]

		#remove genes below the module membership cutoff
		netB_genes = filter(netB_genes, abs(MM) > MM_CUTOFF)$ensembl

		#obtain the ensembl id's of overlapping genes between both study modules
		ensembl_genes = intersect(netA_genes, netB_genes)

		#obtain the HGNC symbols
		hgnc_genes = filter(gene_names, ensembl %in% ensembl_genes)$hgnc

		#obtain the number of overlapping genes
		overlap_len = overlap_mat[i,j]

		#obtain the p-value of overlapping genes
		padj = padj_mat[i,j]

		#if the p-value is statistically signficant
		if (padj < PADJ_CUTOFF){
			#summarize the data as a data frame
			module_overlap = data.frame(ensembl = ensembl_genes,
										hgnc = hgnc_genes,
										module_overlap_len = overlap_len,
										module_overlap_padj = padj)

			#store the results
			sig_overlaps[[netB_module]] = module_overlap
			module_overlaps[[netA_module]][[netB_module]] = module_overlap
		}
	}

	#if list is not empty
	if (length(sig_overlaps) > 0){
		#add a new sheet the excel workbook
		addWorksheet(wb, netA_module)

		#define the starting rows for each data frame in the current sheet
		start_row = 1
		for (k in 2:(length(names(sig_overlaps)) + 1)){
			start_row = c(start_row, ifelse(is.null(nrow(sig_overlaps[[k - 1]])), 0, nrow(sig_overlaps[[k - 1]])) + start_row[k - 1] + 3)
		}

		#define the labels for each data frame
		labels = names(sig_overlaps)

		#write the data frames to the excel worksheet
		for (k in seq_along(sig_overlaps)) {
			writeData(wb, sheet = netA_module, x = labels[k], startRow = start_row[k], startCol = 1)
			writeData(wb, sheet = netA_module, x = sig_overlaps[[k]], startRow = start_row[k] + 1, startCol = 1)
		}
	}
}

#save results as xlsx
saveWorkbook(wb, file = "Significant_module_overlaps.xlsx", overwrite = TRUE)
write(toJSON(module_overlaps), "significant_module_overlaps.json")
