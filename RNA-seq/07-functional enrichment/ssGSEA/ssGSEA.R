#!/usr/bin/Rscript

#Single-sample Gene Set Enrichment Analysis (ssGSEA)
#---
#Script to run an ssGSEA for all samples of an RNA-seq experiment

#define global variables
dir = "/path/to/directory"
norm_read_counts_file = 'your_norm_read_counts.csv' #file name of normalized read counts
biotypes_file = 'your_gene_biotypes.csv' #file name of gene biotypes

count_cutoff = 2
gene_annotation = 'ensembl'
gene_biotype = 'protein_coding' #gene biotypes of interest

#set working directory
setwd(dir)

#load packages
library(dplyr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)


#load data
norm_counts = read.csv(norm_read_counts_file, header = TRUE, sep = ",", row.names = 1) #gene normalized read counts
biotypes = read.csv(biotypes_file) #gene biotypes


#Filter genes by read counts and biotype
#---
#remove genes with average read counts below the count cutoff
norm_counts = norm_counts[rowMeans(norm_counts[]) > count_cutoff,]
norm_counts = data.frame(rownames(norm_counts), norm_counts)
colnames(norm_counts)[1] = gene_annotation

#filter for gene biotypes
norm_counts = merge(biotypes, norm_counts, by = gene_annotation)
norm_counts = filter(norm_counts, biotype == gene_biotype)
norm_counts = column_to_rownames(norm_counts, var = gene_annotation)
norm_counts = norm_counts[,-1]


#Perform ssGSEA for each sample
#---
ssgsea_results = list()

for(sample in colnames(norm_counts)){
	#create the list of genes in decreasing order of normalized read counts
	norm_counts_list = norm_counts[,sample]
	names(norm_counts_list) = rownames(norm_counts)
	norm_counts_list = sort(norm_counts_list, decreasing = TRUE)

	set.seed(123) #set seed to create consistent results (randomness comes from permutations test to asses statistical significance)
	ssgsea = gseGO(geneList = norm_counts_list, #ordered list of genes
                 OrgDb = org.Hs.eg.db, #human DB
                 ont = 'BP', #ontology type
                 keyType = 'ENSEMBL', #gene annotation
                 pAdjustMethod = 'fdr', #adjusted p-value method
				 pvalueCutoff = 0.05, #p-value cutoff
				 minGSSize = 100, #min size of gene set for analyzing
				 maxGSSize = 500, #max size of gene set for analyzing
				 nPermSimple = 1000, #number of permutations for statistical test
				 eps = 0, #boundary for calculating p-values
				 verbose = TRUE)

	ssgsea_results[[sample]] = data.frame(ssgsea)
}

#save results as an excel file
write.xlsx(ssgsea_results, file = "ssGSEA_results.xlsx", rowNames = FALSE)
