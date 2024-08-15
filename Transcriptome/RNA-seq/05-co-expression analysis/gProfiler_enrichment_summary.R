#!/usr/bin/Rscript

#Summarize g:Profiler enrichment results
#---
#Script to summarize the enrichment results obtained from a multiquery output using g:Profiler g:GOSt.

#load packages
library(dplyr)
library(openxlsx)

#define global variables
DIR = "/path/to/directory" #working directory
PADJ_CUTOFF = 0.05

#set working directory
setwd(DIR)

#load data
gProfiler_results = read.csv("gProfiler_raw_multiquery_output.csv")


#Summarize the enrichment results from gProfiler
#---
#define the columns of interest
columns = c("adjusted_p_value__","query_size__","intersection_size__")

#define the module names
modules = colnames(dplyr::select(gProfiler_results, contains(columns)))
modules = unique(sub(".*__", "", modules))

#create an empty list
gProfiler_summary = list()
for (module in modules){
	#subset enrichment results to current module
	gProfiler_result = dplyr::select(gProfiler_results, source, term_id, term_name, term_size, matches(paste0("__", module, "\\b")))
	colnames(gProfiler_result) = c('source','term_id','term_name','term_size','padj','query_size','intersection_size')

	#filter for significant results
	gProfiler_result = filter(gProfiler_result, padj < PADJ_CUTOFF)

	#store results
	gProfiler_summary[[module]] = gProfiler_result
}

#save results as XLSX
write.xlsx(gProfiler_summary, "gProfiler_summary.xlsx")
