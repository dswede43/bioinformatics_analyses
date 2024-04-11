#!/usr/bin/Rscript

#Cytoscape table preparation
#---
#This script prepares the tables needed to create network interaction plots in Cytoscape to visualize
#results from a Gene Set Enrichment Analysis (GSEA). The tables include the matching and name type tables.
#A guide to what these tables are can be found here:
#(https://medium.com/@snippetsbio/how-to-use-cytoscape-for-making-interaction-networks-6-simple-steps-176a1e147020)

#This script also incorperates the log2 Fold Change (LFC) values and p-values produced from differential
#expression analysis.

#define global variables
gsea_df_file = "GSEA_results_table.csv"
LFC_results_file = "LFC_results_table.csv"
gene_names_file = "gene_names_table.csv"
table_outputs = c("matching_table.csv","name_type_table.csv")

#load packages
library(dplyr)
source('https://raw.githubusercontent.com/dswede43/bioinformatics_analyses/main/Functions/create_cytoscape_tables.R')

#load data
gsea_df = read.csv(gsea_df_file) #data frame of GSEA results

LFC_results = read.csv(LFC_results_file) #data frame of LFC results used in the GSEA

gene_names = read.csv(gene_name_file) #data frame of gene annotations


#Create the cytoscape tables
#---
cytoscape_tables = create_cytoscape_tables(gsea_df)


#Matching table
#---
matching_table = cytoscape_tables[[1]]

#convert gene annotations to HGNC symbols
matching_table = merge(matching_table, gene_names[,c(1,3)], by.x = 'gene', by.y = 'ensembl')
matching_table = dplyr::select(matching_table, description, hgnc)
colnames(matching_table) = c('description','gene')


#Name type table
#---
name_type_table = cytoscape_tables[[2]]

#separate name types into its unique types
descriptions_df = filter(name_type_table, type == 'description')
genes_df = filter(name_type_table, type == 'gene')

#append NES values
descriptions_df = merge(descriptions_df, gsea_df[,c('Description','NES')], by.x = 'name', by.y = 'Description')
colnames(descriptions_df) = c('name','type','value')

#append LFC values
genes_df = merge(genes_df, LFC_results, by.x = 'name', by.y = 'ensembl')
colnames(genes_df) = c('name','type','value')

#convert gene annotations to HGNC symbols
genes_df = merge(genes_df, gene_names[,c(1,3)], by.x = 'name', by.y = 'ensembl')
genes_df = dplyr::select(genes_df, hgnc, type, value)
colnames(genes_df) = c('name','type','value')

#merge name types together
name_type_table = rbind(descriptions_df, genes_df)


#Save both tables as csv
#---
write.csv(matching_table, file = table_outputs[1], row.names = FALSE, quote = FALSE)
write.csv(name_type_table, file = table_outputs[2], row.names = FALSE, quote = FALSE)
