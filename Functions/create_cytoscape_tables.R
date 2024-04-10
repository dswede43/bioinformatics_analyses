#Create cytoscape tables
#---
#Function to create the necessary tables needed to create an interaction network in Cytoscape
#(matching and name type tables).
#A guide to what these tables are can be found here:
#(https://medium.com/@snippetsbio/how-to-use-cytoscape-for-making-interaction-networks-6-simple-steps-176a1e147020)

#paramters
#gsea_df: data frame of GSEA results outputted from the gseGO() function from the ClusterProfiler package.

#load packages
library(dplyr)


#Function
#---
create_cytoscape_tables = function(gsea_df){
    
    #create the matchings table
    matchings = data.frame()
    for (i in 1:nrow(gsea_df)){
        
        #obtain the GO term description
        description = gsea_df[i, ]$Description
        
        #obtain the list of associated genes
        genes = strsplit(gsea_df[i, ]$core_enrichment, "/")[[1]]
        
        #create matchings
        matchings = rbind(matchings, data.frame(description = description, gene = genes))
    }
	
	#create the name types table
	name_types = data.frame()
	for (column in colnames(matchings)){
		
		#obtain the unique elements from each column in the matchings data frame
		name_types = rbind(name_types, data.frame(name = unique(matchings[ ,column]), type = column))
	}
	
    return(list(matchings, name_types))
}
