#Append gene names
#---
#Function to append gene annotations to a data frame containing ensembl id's

#paramters
#gene_names: list of gene annotations (HGNC, entrez, GENCODE, etc.) with associated ensembl id's as the list names
#ensembl_df: data frame containing ensembl id's as a column (column must be names 'ensembl')

#load packages
library(dplyr)


#Function
#---
append_gene_names = function(gene_names, ensembl_df){
	
	#convert the gene names to a data frame
	gene_names_df = data.frame(ensembl = names(gene_names), gene_name = gene_names)
	rownames(gene_names_df) = NULL
	
	#join these genes names to the data frame containing ensembl id's
	result = left_join(ensembl_df, gene_names_df, by = 'ensembl')
	
	return(result)
}
