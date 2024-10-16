#!/bin/bash

#SAM file processing - samtools
#---
#This script automates the processing of SAM files including the compression to BAM,
#sorting by genomic coordinates, and indexing using samtools.

#define global variables
DIR="/path/to/directory" #working directory

#create required directories
mkdir -p "${DIR}/outputs/analysis-ready_reads"
mkdir -p "$DIR/outputs/alignments"
mkdir -p "$DIR/logs/alignments"


#Compress outputs to .BAM, sort by genomic coordinates, and index
#---
#create an array of sample names from folder names
mapfile -t sample_names < <(find "${DIR}/samples" -maxdepth 1 -type d ! -path "${DIR}/samples" -exec basename {} \;)

#create empty arrays
for sample_name in ${sample_names[@]}; do
	#record start time
	start_time=$(date +%s)

	#define the input .BAM file
	bam_file="${DIR}/outputs/analysis-ready_reads/${sample_name}.bam"

	#define the sorted output .BAM file
	sorted_bam_file="${DIR}/outputs/analysis-ready_reads/${sample_name}_sorted.bam"

	#sort the .BAM file
	echo "Sorting genomic coordinates for sample ${sample_name}..."
	samtools sort $bam_file -o $sorted_bam_file \
	2> "${DIR}/logs/alignments/samtools_sort-${sample_name}.log"
	echo "Sorting complete!"

	#index the .BAM file
	echo "Indexing for sample ${sample_name}..."
	samtools index $sorted_bam_file \
	2> "${DIR}/logs/alignments/samtools_index-${sample_name}.log"
	echo "Indexing complete!"

	#record end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time))
	minutes=$(( run_time / 60 ))
    echo "Sample ${sample_name} took ${minutes} minutes"
done
