#####RNA-SEQ SIMULATION
#The following script simulates the mapped read counts from an RNA-sequencing experiment using a negative binomial distribution,
#where the mean and dispersion estimates are taken from a real life RNA-seq experiment.
#Differentially expressed (DE) genes are manually created to test the specificity and sensitivity of various DE methods.

#load packages
library(DESeq2)
library(matrixStats)
library(ggplot2)
library(ggpubr)

#define file directory
setwd("C:/path/to/data")

#define global variables
data_file = "mapped_read_counts.csv" #define read count matrix file name
n_genes = 20000 #number of genes
n_samples = 10 #sample size


####GENERATE RNA-SEQ DATA CHARACTERISTICS (MEAN AND VARIANCE) FROM REAL DATA
#load RNA-seq data
real_data = read.csv(data_file, header = TRUE, sep = ",", row.names = 1) #read count matrix from real data
real_data = real_data[(rowSums(real_data != 0) > 0),] #remove all rows with zeroes across all columns


###CREATE RANDOM SUBSET OF GENES
set.seed(123) #set seed
gene_subset = sample(nrow(real_data), size = n_genes, replace = FALSE) #create subset of random genes from real data
real_subset = real_data[gene_subset,] #subset read count data


###ESTIMATE PARAMETERS FOR EACH GENE (MEAN, VARIANCE, AND DISPERSION)
gene_means = rowMeans(real_subset) #estimate read count mean for each gene
gene_vars = rowVars(as.matrix(real_subset), useNames = TRUE) #estimate read count variance for each gene

#estimate read count dispersion for each gene using DESeq2
meta = data.frame(sample = colnames(real_data)) #create a metadata table with sample names (required for DESeq2)
rownames(meta) = colnames(real_data) #set rownames of metadata
dds = DESeqDataSetFromMatrix(countData = real_subset, #create a DESeq2 object
							 colData = meta,
							 design = ~ 1)
dds = estimateSizeFactors(dds) #estimate gene size factors
dds = estimateDispersions(dds) #estimate gene dispersions
gene_disps = dispersions(dds)

#create table of gene parameters
real_params = data.frame(mean = gene_means, var = gene_vars, disp = gene_disps)
####


####GENERATE SIMULATED RNA-SEQ DATA
#define a function to simulate counts for a single gene
simulate_counts <- function(mean, disp, n_samples){
	return(rnbinom(n = n_samples, mu = mean, size = disp)) #simulate read counts using a negative binomial distribution
}


###CONTROL SAMPLE SIMULATION
#use lapply to simulate counts for each gene
tmp = lapply(1:nrow(real_params), function(i){
	simulate_counts(real_params[i, 1], real_params[i, 3], n_samples)
})
control_data = as.data.frame(do.call(rbind, tmp)) #convert the list to a data frame


###TREATMENT SAMPLE SIMULATION
#create differentially expressed (DE) gene for the treatment sample simulation
#~10% of genes are DE
	#5% up-regulated
		#3% 1.5x
		#1.5% 2x
		#0.5% 3x
	#5% down-regulated
		#3% 1.5x
		#1.5% 2x
		#0.5% 3x

set.seed(123) #set seed
tmp = sample(nrow(real_params), size = round(0.1 * nrow(real_params)), replace = FALSE) #take a 10% random sample of genes
DE_params = real_params[-tmp,] #remove those sampled genes from the gene parameter data frame
DE_genes = real_params[tmp,] #create a new data frame with the parameters of those sampled genes

#create DE genes through multiplying gene means
DE_genes[1:600, "mean"] = DE_genes[1:600, "mean"] * 3/2 #1.5x up-regulated
DE_genes[601:1200, "mean"] = DE_genes[601:1200, "mean"] * 2/3 #1.5x down-regulated
DE_genes[1201:1500, "mean"] = DE_genes[1201:1500, "mean"] * 2 #2x up-regulated
DE_genes[1501:1800, "mean"] = DE_genes[1501:1800, "mean"] * 1/2 #2x down-regulated
DE_genes[1801:1900, "mean"] = DE_genes[1801:1900, "mean"] * 3 #3x up-regulated
DE_genes[1901:2000, "mean"] = DE_genes[1901:2000, "mean"] * 1/3 #3x down-regulated
DE_params = rbind(DE_genes, DE_params) #row bind DE gene parameters to original gene parameters

#use lapply to simulate counts for each gene
tmp = lapply(1:nrow(DE_params), function(i){
	simulate_counts(DE_params[i, 1], DE_params[i, 3], n_samples)
})
treatment_data = as.data.frame(do.call(rbind, tmp)) #convert the list to a data frame

sim_data = data.frame(control_data, treatment_data) #create the simulated data frame

#add row and column names to simulated data
control_names = paste0("control_", 1:n_samples) #create control sample names
treatment_names = paste0("treatment", 1:n_samples) #create treatment sample names
gene_names = paste0("gene", 1:n_genes) #create gene names
colnames(sim_data) = c(control_names, treatment_names)
rownames(sim_data) = gene_names

#save simulated data as a csv
write.csv(sim_data, "Simulated_RNA-seq_data.csv")
####


####VISUALIZE ALL DATA
###VISUALIZE THE REAL DATA
#gene read count distribution
real_data_hist = ggplot(data = real_params, aes(x = mean)) +
	geom_histogram(bins = 50) +
	labs(x = "mean read counts",
		 title = paste0("Real data: histogram of mean\nread counts for ", n_genes, " genes")) +
	theme_classic()

#gene read count mean-variance scatterplot
real_data_plot = ggplot(data = real_params, aes(x = log10(mean), y = log10(var))) +
	geom_point(size = 1, alpha = 0.05) +
	geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
	labs(x = "log10(mean read counts)", y = "log10(variance of read counts)",
		 title = paste0("Real data: mean-variance relationship\nof read counts for ", n_genes, " genes")) +
	theme_classic()


###VISUALIZE THE SIMULATED DATA
#calculate mean and variance of simulated data
sim_params = data.frame(mean = rowMeans(sim_data),
						var = rowVars(as.matrix(sim_data), useNames = TRUE))

#gene read count distribution
sim_data_hist = ggplot(data = real_params, aes(x = mean)) +
	geom_histogram(bins = 50) +
	labs(x = "mean read counts",
		 title = paste0("Simulated data: histogram of mean\nread counts for ", n_genes, "genes")) +
	theme_classic()

#gene read count mean-variance scatterplot
sim_data_plot = ggplot(data = sim_params, aes(x = log10(mean), y = log10(var))) +
	geom_point(size = 1, alpha = 0.05) +
	geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
	labs(x = "log10(mean read counts)", y = "log10(variance of read counts)",
		 title = paste0("Simulated data: mean-variance relationship\nof read counts for ", n_genes, " genes")) +
	theme_classic()

#combine all visuals into one plot
plots = ggarrange(real_data_hist, real_data_plot,
				  sim_data_hist, sim_data_plot,
				  nrow = 2, ncol = 2)

#save plot as pdf
pdf("Real_&_simulated_data_plots.pdf", height = 6, width = 8)
plots
dev.off()
####
#####