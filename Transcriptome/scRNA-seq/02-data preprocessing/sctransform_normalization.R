##!/usr/bin/Rscript

#Seurat sctransform
#---
#This script utilizes the sctransform framework from Seurat to normalize and transform scRNA-seq read count data.
#This framework replaces the NormalizeData(), ScaleData(), and FindVariableFeatures() functions typically done
#in a standard scRNA-seq preprocessing workflow.

#load packages
library(ggplot2)
library(Seurat)
library(sctransform)

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


#Apply sctransform normalization
#---
#calculate the percentage of counts originating from mitochondrial genes
seurat_obj = PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")

#run sctransform
seurat_obj = SCTransform(seurat_obj, #specify seurat object
						 vars.to.regress = "percent.mt", #control for variability associated with mitochondrial contamination
						 verbose = TRUE)


#Save seurat object as a .RDS
#---
saveRDS(seurat_obj, file = seurat_output_file)
