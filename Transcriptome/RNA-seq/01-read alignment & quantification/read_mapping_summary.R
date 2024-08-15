#!/usr/bin/Rscript

#Read alignment and quantification summary
#---
#This script summarizes and visualizes the count results of read alignment and quantification.

#required directory structure:
#---mappings/
#---quants/

#load packages
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)

#define global variables
dir = "/path/to/directory" #working directory

#set working directory
setwd(dir)

#define the sample names
samples = list.dirs(path = "quants/", full.names = FALSE)
samples = samples[samples != ""]


#Alignment results
#---
alignments = list()
for (sample in samples){
    #load the log alignment file
    log_file = read.csv(paste0("mappings/", sample, "/alignment.log"))
    colnames(log_file) = "counts"

    #extract count statistics
    counts = c()
    for(i in 1:4){
        #extract the numbers larger than 1-digit and with white space to the right
        counts = c(counts, as.numeric(regmatches(log_file[i,], gregexpr("\\b\\d{2,}\\b+\\s", log_file[i,]))))
    }

    #store the results
    alignments[[paste0("S", sample)]] = counts
}

#convert results into a data frame
alignments = data.frame(t(data.frame(alignments)))
colnames(alignments) = c("total_reads","unaligned","aligned_once","aligned_multiple")

#save results as CSV
write.csv(alignments, "mappings/alignment_stats.csv", row.names = TRUE)


#Quantification results
#---
quants = list()
for(sample in samples){
    #load the quantification summary file
    summary_file = read.csv(paste0("quants/", sample, "/sample", sample, ".txt.summary"), sep = "\t", row.names = 1)
    colnames(summary_file) = "counts"

    #extract count statistics
    counts = summary_file[c("Assigned","Unassigned_Unmapped","Unassigned_MultiMapping","Unassigned_NoFeatures","Unassigned_Ambiguity"),]

    #store the results    
    quants[[paste0("S", sample)]] = counts
}

#convert results into a data frame
quants = data.frame(t(data.frame(quants)))
colnames(quants) = c("Assigned","Unassigned_Unmapped","Unassigned_MultiMapping","Unassigned_NoFeatures","Unassigned_Ambiguity")

#save results as CSV
write.csv(quants, "quants/quantification_stats.csv", row.names = TRUE)


#Visualizations
#---
#visualizations of alignment and quantification results


#Alignment results
#---
#pivot the results into long format
alignments_df = rownames_to_column(alignments, var = "sample")
alignments_df = alignments_df[-2]
alignments_df = pivot_longer(alignments_df, cols = unaligned:aligned_multiple, names_to = "alignment", values_to = "count")

#calculate the proportions
alignments_df = alignments_df %>%
    mutate(sample = factor(sample), alignment = factor(alignment)) %>%
    group_by(sample) %>%
    mutate(prop = count / sum(count))

#stacked bar plot of proportions
align_prop = ggplot(data = alignments_df, aes(x = sample, y = prop, fill = alignment)) +
	geom_col(aes(fill = alignment), position = position_stack(), col = "black") +
	geom_text(aes(label = round(prop, 2)), size = 2, position = position_stack(vjust = 0.5)) +
	labs(x = "Sample", y = "Proportions") +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_discrete(limits = rev) +
	coord_flip() +
	theme_classic()

#stacked bar plot of counts
align_count = ggplot(data = alignments_df, aes(x = sample, y = count, fill = alignment)) +
	geom_col(aes(fill = alignment), position = position_stack(), col = "black") +
	geom_text(aes(label = round(count / 10^{6}, 2)), size = 2, position = position_stack(vjust = 0.5)) +
	labs(x = "Sample", y = "Counts (millions)") +
	scale_y_continuous(expand = c(0,0),
		limits = c(0,100000000),
		breaks = c(0,25000000,50000000,75000000,100000000),
		labels = c("0","25","50","75","100")) +
	scale_x_discrete(limits = rev) +
	coord_flip() +
	theme_classic()


#Quantification results
#---
#pivot the results into long format
quants_df = rownames_to_column(quants, var = "sample")
quants_df = dplyr::select(quants_df, sample, Assigned, Unassigned_NoFeatures, Unassigned_Ambiguity) #sum of these columns add up to the number of reads aligned once
quants_df = pivot_longer(quants_df, cols = Assigned:Unassigned_Ambiguity, names_to = "quantification", values_to = "count")

#calculate the proportions
quants_df = quants_df %>%
	mutate(sample = factor(sample), quantification = factor(quantification)) %>%
	group_by(sample) %>%
	mutate(prop = count/sum(count))

#stacked bar plot of proportions
quant_prop = ggplot(data = quants_df, aes(x = sample, y = prop, fill = quantification)) +
	geom_col(aes(fill = quantification), position = position_stack(), col = "black") +
	geom_text(aes(label = round(prop, 2)), size = 2, position = position_stack(vjust = 0.5)) +
	labs(x = "Sample", y = "Proportions") +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_discrete(limits = rev) +
	coord_flip() +
	theme_classic()

#stacked bar plot of counts
quant_count = ggplot(data = quants_df, aes(x = sample, y = count, fill = quantification)) +
	geom_col(aes(fill = quantification), position = position_stack(), col = "black") +
	geom_text(aes(label = round(count / 10^{6}, 2)), size = 2, position = position_stack(vjust = 0.5)) +
	labs(x = "Sample", y = "Counts (millions)") +
	scale_y_continuous(expand = c(0,0),
		breaks = c(0,10000000,20000000,30000000,40000000,50000000,60000000),
		labels = c("0","10","20","30","40","50","60")) +
	scale_x_discrete(limits = rev) +
	coord_flip() +
	theme_classic()

#save the plots as pdf
plots = list(align_prop, align_count, quant_prop, quant_count)
pdf("mapping_summary_plots.pdf", height = 8, width = 8)
ggarrange(plotlist = plots, ncol = 2, nrow = 2)
dev.off()
