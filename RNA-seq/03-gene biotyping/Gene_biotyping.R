#!/usr/bin/Rscript

#Gene naming and biotyping
#---
#Script to associate a different gene names and biotypes for a list of ensembl genes using the Ensembl database via Biomart.

#define global variables
dir = "/path/to/directory"
read_count_file = 'your_read_counts.csv' #read count matrix file name
organism_db = 'hsapiens_gene_ensembl' #define the organism to query

#set working directory
setwd(dir)

#load packages
library(biomaRt)

#load data
#gene ensembl id's
ensembl_genes = rownames(read.csv(read_count_file, header = TRUE, sep = ",", row.names = 1))


#Gene biotyping
#---
#connect to the Ensembl database via biomart
mart = useMart("ensembl", dataset = organism_db)

#search the available attributes in this Mart
searchAttributes(mart = mart, pattern = 'hgnc')
searchAttributes(mart = mart, pattern = 'entrez')
searchAttributes(mart = mart, pattern = 'biotype')

#define the attributes to query
attributes = c('hgnc_symbol','entrezgene_id','gene_biotype')

#query for gene biotypes
gene_names_biotypes = getBM(filters = "ensembl_gene_id", #define qeury input
                      attributes = attributes, #define qeury output
					  values = ensembl_genes, #list of ensembl genes
					  mart = mart)
colnames(gene_names_biotypes) = c('ensembl','hgnc','entrez','biotype')

#save data frame as csv
write.csv(gene_names_biotypes, "gene_names_biotypes.csv", row.names = TRUE)
