#Filter gene biotype
#---
#Function to filter a data frame containing gene ensembl id's to a specific gene biotype

#paramters
#gene_biotypes: list of gene biotypes (protein-coding, lncRNA, snRNA, etc.) with associated ensembl id's as the list names
#ensembl_df: data frame containing ensembl id's as a column (column must be named 'ensembl')
#biotype: name of the biotype to keep in the 'ensembl_df' data frame


#load packages
library(dplyr)


#Function
#---
filter_gene_biotype = function(gene_biotypes, ensembl_df, biotype = 'protein_coding'){
	
	#convert the gene names to a data frame
	gene_biotypes_df = data.frame(ensembl = names(gene_biotypes), gene_biotype = gene_biotypes)
	rownames(gene_biotypes_df) = NULL
	
	#join these genes names to the data frame containing ensembl id's
	result = left_join(ensembl_df, gene_biotypes_df, by = 'ensembl')
	
	#remove non-protein coding biotypes
	result = filter(result, gene_biotype == biotype)
	
	return(result)
}
