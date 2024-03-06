#Semantic similarity reduction
#---
#Function to remove redundant GO terms based on semantic similarity and normalized enrichment scores (NES)

#parameters
#GO_enr: list of normalized enrichment scores (NES) calculated from GSEA with their GO term id names
#C: semantic similarity cutoff

#load packages
library(GOSim)
library(dplyr)


#Function
#---
#define a function to remove redundant GO terms
reduce_GO_list = function(GO_enr, C){
	
	#define list of GO terms to reduce
	GO_terms = names(GO_enr)
	
	#calculate the pairwise similarity matrix between GO terms
	sim_mat = getTermSim(GO_terms, method = 'Lin')
	
	#ensure the max semantic similarity value is 1
	sim_mat[sim_mat > 1] = 1
	
	GO_comparisons = data.frame()
	#for each value in the lower semantic similarity matrix
	for (i in 1:(nrow(sim_mat) - 1)){
		for (j in (i + 1):ncol(sim_mat)){
			
			#if the semantic similarity score is greater than the cutoff
			if (sim_mat[i,j] > C){
				
				#obtain the GO term id's for this comparison
				GO_comparisons = rbind(GO_comparisons, c(rownames(sim_mat)[i], colnames(sim_mat)[j]))
			}
		}
	}
	
	#if comparisons greater than the cutoff exist
	if (nrow(GO_comparisons) != 0){
		
		#add column headers to GO comparisons
		colnames(GO_comparisons) = c('GO1','GO2')
	
		GO_remove = c()
		#for each GO term comparison
		for (i in 1:nrow(GO_comparisons)){
		
			#obtain the GO term enrichment scores
			GO_comparison = as.character(GO_comparisons[i,])
			GO_enr_comparison = GO_enr[names(GO_enr) %in% GO_comparison]

			#compare the absolute value of enrichment scores
			if (abs(GO_enr_comparison[1]) < abs(GO_enr_comparison[2])){
				
				#retain the GO id with the lesser enrichment score
				GO_remove = c(GO_remove, names(GO_enr_comparison[1]))
			} else{
				GO_remove = c(GO_remove, names(GO_enr_comparison[2]))
			}
		}
	}
	else{
		GO_remove = NULL
	}
	
	#remove the redundant GO terms from the original list
	GO_reduced = rownames(sim_mat)[!rownames(sim_mat) %in% unique(GO_remove)]
	
	#return the list of reduced GO terms
	return(GO_reduced)
}
