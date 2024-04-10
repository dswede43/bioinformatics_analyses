#Subset GO terms
#---
#Function to subset GO terms to those associated with a list of specific genes.

#paramters
#gsea_df: data frame of GSEA results outputted from the gseGO() function from the ClusterProfiler package.
#ensembl_list: list of ensembl annotated genes.


#Function
#---
subset_GO_terms = function(gsea_df, ensembl_list){
	GO_ids = c()
	for (i in 1:nrow(gsea_df)){
		
		#create the list of genes associated to the GO term
		ensembl_GO = strsplit(gsea_df[i, ]$core_enrichment, "/")[[1]]
		
		#if there is a non-null intersection between gene lists
		if(length(intersect(ensembl_GO, ensembl_list)) > 0){
			
			#append the current GO term to the list
			GO_ids = c(GO_ids, gsea_df[i, ]$ID)
		}
	}
	
	#filter the GSEA data frame by the subset of GO terms
	gsea_subset = filter(gsea_df, ID %in% GO_ids)
	
	return(gsea_subset)
}
