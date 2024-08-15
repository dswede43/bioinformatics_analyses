##!/usr/bin/Rscript

#Linear dimension reduction (principal component analysis)
#---
#This script completes a linear dimension reduction of scRNA-seq scaled count data using principal component analysis.

#load packages
library(Seurat)
library(ggplot2)

#define global variables
sparse_mat_dir = "PBMC_data/"
project_name = "PBMC"
seurat_output_file = paste0(sparse_mat_dir, "PBMC_seurat_preprocessed.rds")


#Load the dataset
#---
#set the working directory of data files
sparse_mat = Read10X(data.dir = sparse_mat_dir)

#create the Seurat object
seurat_obj = CreateSeuratObject(counts = sparse_mat, #unnormalized counts
								project = project_name, #name of project
								min.cells = 3, #include features detected in at least this many cells
								min.features = 200) #include cells detected in at least this many features


#Perform linear dimensional reduction
#---
seurat_obj = RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

#visualize the loading scores from PCA
VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")

#visualize the first two principal components
DimPlot(seurat_obj, reduction = "pca") + NoLegend()

#visualize sources of gene heterogeneity across cells
DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)


#Determine the number of principal components
#---
#elbow plot method (PC standard deviations)
ElbowPlot(seurat_obj)

#elbow bar graph method (proportional variance explained by PCs)
#calculate proportion of variance explained by PCs
var = seurat_obj[["pca"]]@stdev^2
per = round(var / sum(var)*100, 1)

#create the data frame
pca_per = data.frame(per = per[1:20], PC = seq(1:20))

#plot the bar graph
ggplot(data = pca_per, aes(x = PC, y = per)) +
	geom_bar(stat = "identity") +
	labs(x = "principal component", y = "proportion of variance explained (%)") +
	theme_classic()

#define the number of PCs
pc_num = 5


#Cell clustering
#---
#contruct the K-Nearest Neigbour (KNN) graph
seurat_obj = FindNeighbors(seurat_obj, dims = 1:pc_num)

#cluster the cells
seurat_obj = FindClusters(seurat_obj, resolution = 0.5)


#Principal component plot
#---
#create the data frame of PC scores and cell clusters
pca_scores = data.frame(seurat_obj[["pca"]]@cell.embeddings, cluster = Idents(seurat_obj))

#visualize the first two PCs
ggplot(data = pca_scores, aes(x = PC_1, y = PC_2, color = cluster)) +
	geom_point(size = 1) +
	labs(x = paste0('PC1 - ', per[1], '%'),
		 y = paste0('PC2 - ', per[2], '%')) +
	theme_classic()


#Save seurat object as a .RDS
#---
saveRDS(seurat_obj, file = seurat_output_file)
