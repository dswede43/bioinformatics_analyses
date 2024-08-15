#!/usr/bin/Rscript

#Query differential expression (DE) results from DESeq2
#---
#Script to query the results table from a differential expression (DE) analysis with DESeq2.

#Table of contents:
#---
#1. Query the differential expressed genes
#2. Query a subset of differential expressed genes

#load packages
library(dplyr)
library(jsonlite)
library(openxlsx)

#define global variables
DIR = "/path/to/directory" #working directory
LFC_CUTOFF = 0.5 #log2 fold change cutoff
PADJ_CUTOFF = 0.1 #adjusted p-value cutoff

#set working directory
setwd(DIR)

#load data
DE_results_table = fromJSON("DESeq2_results_table.JSON")
gene_biotypes = read.csv("RNA_biotypes.csv")
gene_names = read.csv("gene_annotations.csv")
gene_subset = read.delim("gene_subset.txt", header = FALSE)$V1


#1. Query the differential expressed genes
#---
#filter for DE genes
DE_result = filter(DE_results_table, abs(log2FoldChange) > LFC_CUTOFF & padj < PADJ_CUTOFF)

#merge HGNC symbols
DE_result = merge(DE_result, gene_names, by = 'ensembl')

#merge gene biotypes
DE_result = merge(DE_result, gene_biotypes, by = 'ensembl')

#store the results
DE_result = dplyr::select(DE_result, ensembl, hgnc, biotype, log2FoldChange, padj)

#save results as CSV
write.csv(DE_result, "DE_gene_results.csv")


#2. Query a subset of differential expressed genes
#---
DE_subset = filter(DE_result, ensembl %in% gene_subset)

#save results as CSV
write.csv(DE_subset, "DE_gene_subset_results.csv")
