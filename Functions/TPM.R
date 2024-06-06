#Transcripts Per Million (TPM) normalization
#---
#Function to normalize RNA-seq read counts using Transcripts Per Million (TPM).

#parameters
#read_counts: matrix of read counts where rows are genes and columns are samples
#gene_lengths: max lengths of each gene


#Function
#---
tpm = function(read_counts, gene_lengths){
	#normalize for gene length (within-sample effect)
	x = read_counts / gene_lengths

	#normalize for read-depth (between-sample effect)
	return (t(t(x) * 1e6 / colSums(x)))
}
