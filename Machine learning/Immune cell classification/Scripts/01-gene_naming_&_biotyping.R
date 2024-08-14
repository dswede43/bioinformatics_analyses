#!/usr/bin/Rscript

#Gene naming and biotyping
#---
#Script to associate a different gene names and biotypes for a list of ensembl genes using the Ensembl database via Biomart.


#load packages
library(biomaRt)

#define global variables
DIR = "C:/Users/dstra/Nextcloud/Documents/Computers/Code/Github repos/bioinformatics analyses/Machine learning/Immune cell classification" #working directory
ORGANISM_DB = 'hsapiens_gene_ensembl' #define the organism to query

#set working directory
setwd(DIR)

#load data
ensembl_genes = read.csv("Data/TPM_gene_counts.csv")$ensembl #gene ensembl id's


#Gene biotyping
#---
#connect to the Ensembl database via biomart
mart = useMart("ensembl", dataset = ORGANISM_DB)

#define the attributes to query
attributes = c('hgnc_symbol','gene_biotype')

#query for gene biotypes
gene_names_biotypes = getBM(filters = "ensembl_gene_id", #define qeury input
							attributes = attributes, #define qeury output
							values = ensembl_genes, #list of ensembl genes
							mart = mart) #biomart database
colnames(gene_names_biotypes) = c('ensembl','hgnc','biotype')

#save data frame as csv
write.csv(gene_names_biotypes, "Data/gene_names_biotypes.csv", row.names = TRUE)
