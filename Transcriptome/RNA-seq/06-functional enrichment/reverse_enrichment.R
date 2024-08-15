#!/usr/bin/Rscript

#Reverse enrichment analysis
#---
#Script to complete a reverse enrichment analysis. ie. filter for specific GO terms that map to specific gene sets
#and contain certain keywords in their description.


#Table of contents:
#---
#1. Obtain the GO terms
#2. Plot GO heatmaps
#3. Plot gene heatmaps


#load packages
library(dplyr)
library(jsonlite)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)

#define global variables
DIR = "/path/to/directory" #working directory
COMPARISONS = c("BDCvsHDT2","HDT2vsHDT30","HDT2vsHDT60","HDT60vsR1","R1vsR12","R1vsR30")
KEYWORDS = c("virus","viral") #keywords to query
CUTOFFS = c(0.5, 0.05)

#set working directory
setwd(DIR)

#load data
DE_results = fromJSON("DESeq2_results.json")
gene_names = read.csv("gene_identities.csv")[ , c(1,3)]


#Functions
#---
#define function to return GO term profiles for a gene set
geneset_GO_profile = function(geneset){
	#Function to return GO term profiles for a geneset across all ontology levels.

	#Arguments:
	#geneset -- array of gene identities (such as ensembl).

	#Returns:
	#GO_profiles -- data frame of GO profile results from the gene set.

	#define the GO ontologies and their max levels
	onts = list("BP" = 19,"MF" = 14,"CC" = 15)

	#define an empty data frame
	GO_profiles = data.frame()
	for (i in 1:length(onts)){
		#define the current ontology and it's max level
		ont = names(onts[i])
		max_level = onts[[i]]

		for (level in 1:max_level){
			#obtain the GO profile
			GO_profile = groupGO(
				geneset,
				OrgDb = org.Hs.eg.db,
				ont = ont,
				keyType = "ENSEMBL",
				level = level,
				readable = FALSE)

			#convert to data frame
			GO_profile = data.frame(GO_profile)

			#append current ontology
			GO_profile$ontology = ont

			#store the results
			GO_profiles = rbind(GO_profiles, GO_profile)
		}
	}
	#remove empty gene mappings
	GO_profiles = GO_profiles[!GO_profiles$geneID == "", ]

	#remove duplicate terms
	GO_profiles = unique(GO_profiles)

	return(GO_profiles)
}

#define function to calculate the median LFCs for each GO terms geneset
get_gene_LFC_medians = function(geneIDs, LFCs){
	#Function to calculate the median log2FoldChanges for each GO terms geneset.

	#Arguments:
	#geneIDs -- array of '/' seperated gene identities obtained after mapping them
		#to their associated GO terms using clusterProfiler groupGO() function.
	#LFCs -- array of gene LFC values along with the gene identity names.

	#Returns:
	#LFC_medians -- array of calculated LFC medians.

	#define an empty array
	LFC_medians = c()
	for (i in 1:length(geneIDs)){
		#define the current GO terms geneset
		geneset = geneIDs[i]
		geneset = strsplit(geneset, split = "/")[[1]]

		#calculate the median LFC for this geneset
		LFC_median = median(LFCs[geneset])

		#store the results
		LFC_medians = c(LFC_medians, LFC_median)
		}
	return(LFC_medians)
}


#1. Obtain the GO terms
#---
#define an empty array
geneset = c()
for (comparison in COMPARISONS){
	#define the DE results
	DE_comparison = DE_result[[comparison]]

	#define the DE gene identities
	DE_genes = filter(DE_comparison, abs(log2FoldChange) > cutoffs[1] & padj < cutoffs[2])$ensembl

	#store the results
	geneset = c(geneset, DE_genes)
}

#obtain the unique geneset
geneset = unique(geneset)

#obtain the DE geneset GO profile
GO_profiles = geneset_GO_profile(geneset)

#filter the GO profile for specific terms
GO_filtered = filter(GO_profiles, grepl(paste(KEYWORDS, collapse = '|'), Description))

#define the GO term geneset identities
geneIDs = GO_filtered$geneID

for (comparison in COMPARISONS){
	#define the DE results and LFC values
	DE_comparison = DE_result[[comparison]]
	LFCs = DE_comparison$log2FoldChange
	names(LFCs) = DE_comparison$ensembl

	#calculate the median LFC for each GO term geneset
	LFC_medians = get_gene_LFC_medians(geneIDs, LFCs)

	#store the results
	GO_filtered[ , comparison] = LFC_medians
}

#save results as CSV
write.csv(GO_filtered, "enriched_GO_terms.csv", row.names = FALSE)


#2. Plot GO heatmaps
#---
#define the LFC cutoff
LFC_cutoff = 0.5

#define the heatmap aesthetics
breaks = seq(-1.5, 1.5, by = 0.05)
colors = colorRampPalette(c("#FFDD00","#FFFFFF","#00BBFF"))(length(breaks))

#setup the plotting data
plot_data = GO_filtered[ , c(2,6:12)]
plot_data$label = paste0("GO:", plot_data$ontology, " - ", plot_data$Description)
rownames(plot_data) = plot_data$label
plot_data = plot_data[ , -c(1:2,9)]

#remove terms below the LFC cutoff
plot_mat = plot_data[rowSums(plot_data > LFC_cutoff) >= 1, ]

#create the heatmap
pheatmap(
	plot_mat,
	border_cols = FALSE,
	cluster_cols = FALSE,
	cluster_rows = TRUE,
	angle_col = 45,
	display_numbers = TRUE,
	color = colors,
	breaks = breaks,
	cellheight = 12,
	cellwidth = 20,
	fontsize = 8,
	height = 6,
	width = 7,
	filename = paste0("reverse_enrichment_GO_term_heatmap.pdf"))


#3. Plot gene heatmaps
#---
#obtain the GO terms above the LFC cutoff
tmp = GO_filtered[ , c(7:12)]
terms = rownames(tmp[rowSums(tmp > LFC_cutoff) >= 1, ])
plot_data = filter(GO_filtered, ID %in% terms)

#find all gene identities associated to each GO term
geneset = c()
for (i in 1:nrow(plot_data)){
	#define the current GO terms geneset
	genes = plot_data$geneID[i]
	genes = strsplit(genes, split = "/")[[1]]

	#store the results
	geneset = c(geneset, genes)
}

#remove duplicate genes
geneset = unique(geneset)

#define the current studies DE results
DE_result = DE_results[[study]]

#define the current studies DE comparisons
comparisons = COMPARISONS[[study]]

#obtain the gene LFC values
plot_list = list()
for (comparison in comparisons){
	#define the DE results
	DE_comparison = DE_result[[comparison]]

	#obtain the LFC values
	LFCs = filter(DE_comparison, ensembl %in% geneset)[ , c("ensembl","log2FoldChange")]

	#store the results
	plot_list[[comparison]] = LFCs

#setup the plotting matrix
plot_mat = do.call(cbind, plot_list)
colnames(plot_mat)[1] = "ensembl"
plot_mat = merge(plot_mat, gene_names, by = 'ensembl')
rownames(plot_mat) = plot_mat$hgnc
plot_mat = plot_mat[ , grep("log2FoldChange", colnames(plot_mat))]
colnames(plot_mat) = comparisons

#create the heatmap
pheatmap(
	plot_mat,
	border_cols = FALSE,
	cluster_cols = FALSE,
	cluster_rows = TRUE,
	angle_col = 45,
	display_numbers = TRUE,
	color = colors,
	breaks = breaks,
	cellheight = 12,
	cellwidth = 20,
	fontsize = 8,
	height = 7,
	width = 5,
	filename = paste0("reverse_enrichment_gene_heatmap.pdf"))
