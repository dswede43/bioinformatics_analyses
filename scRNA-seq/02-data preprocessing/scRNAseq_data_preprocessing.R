##!/usr/bin/Rscript

#Seurat data preprocessing
#---
#This script completes the standard preprocessing workflow for a single-cell RNA-sequencing experiment.
#This includes cell filtration, data normalization/scaling, and detection of highly variable genes (feature selection).

#load packages
library(dplyr)
library(Seurat)

#define global variables
unique_count_thresholds = c(200, 2500) #number of unique gene counts thresholds
mt_count_perc = 5 #threshold for the number of reads originating from mitochondrial genes
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


#Quality control of cells
#---
#calculate the percentage of counts originating from mitochondrial genes
#(low quality/dying cells exhibit mitochondrial contamination)
seurat_obj[["percent.mt"]] = PercentageFeatureSet(seurat_obj, pattern = "^MT-")

#filter cells based on unique gene count and mitochondrial percentage thresholds
seurat_obj = subset(seurat_obj,
					subset = nFeature_RNA > unique_count_thresholds[1] &
							 nFeature_RNA < unique_count_thresholds[2] &
							 percent.mt < mt_count_perc)


#Read normalization
#---
seurat_obj = NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)


#Identification of highly variable genes (feature selection)
#---
#use the variance stablization transformation (VST) to identify highly variable genes
seurat_obj = FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

#obtain the top 10 most variable genes
top10 = head(VariableFeatures(seurat_obj), 10)

#plot the results
plot = VariableFeaturePlot(seurat_obj)
plot = LabelPoints(plot = plot, points = top10, repel = TRUE)


#Scaling data
#---
#define the genes to scale
all_genes = rownames(seurat_obj)

#apply scaling
seurat_obj = ScaleData(seurat_obj, features = all_genes)


#Save seurat object as a .RDS
#---
saveRDS(seurat_obj, file = seurat_output_file)
