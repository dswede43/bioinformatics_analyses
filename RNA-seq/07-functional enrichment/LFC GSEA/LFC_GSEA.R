#!/usr/bin/Rscript

#LFC GSEA
#---
#This script completes Gene Set Enrichment Analyses (GSEA) using the log2 fold changes (LFC)
#obtained from each differential expression analysis comparison. Each GSEA is handled by the
#clusterProfiler package (https://guangchuangyu.github.io/software/clusterProfiler/).

#define global variables
dir = "/path/to/directory" #directory to save data
LFC_results_file = 'your_LFC_results.csv' #file name of LFC results from differential expression analysis
biotypes_file = 'your_gene_biotypes.csv' #file name of gene biotypes

gene_header = 'ensembl' #gene annotation name
gene_biotype = 'protein_coding' #gene biotypes of interest
comparisons = c('your','differential','expression','comparisons') #define differential expression (DE) comparisons (***must be the same as column header names***)

#set working directory
setwd(dir)

#load packages
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)

#load data
LFC_results = read.csv(LFC_results_file, header = TRUE, sep = ",") #LFC results from differential expression analysis
biotypes = read.csv(biotypes_file, header = TRUE, sep = ",") #ensembl RNA biotypes


#subset data to specific gene biotypes
#---
LFC_df = LFC_results[,c(gene_header, comparisons)]
LFC_df = merge(biotypes, LFC_df, by = 'ensembl')
LFC_df = filter(LFC_df, biotype == gene_biotype)


#Perform GSEA of LFC's
#---
gsea_results = list()

for(comparison in comparisons){
	#create the list of genes in decreasing order of LFC
    LFC_list = LFC_df[,comparison]
    names(LFC_list) = LFC_df[,gene_header]
    LFC_list = sort(LFC_list, decreasing = TRUE)

	set.seed(123) #set seed to create consistent results (randomness comes from permutations test to asses statistical significance)
    gsea = gseGO(geneList = LFC_list, #ordered list of genes
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

	gsea_results[[comparison]] = data.frame(gsea)
}

#save an example GSEA plot to illustrate how enrichment scores are calculated
pdf("GSEA_enrichment_score_example.pdf", width = 8, height = 8)
gseaplot(gsea, gsea_results[[comparison]][1,1])
dev.off()

#save results as an excel file
write.xlsx(gsea_results, file = "LFC_GSEA_results.xlsx", rowNames = F)
