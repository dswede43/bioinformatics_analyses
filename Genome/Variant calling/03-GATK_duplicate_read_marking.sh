#!/bin/bash

#Marking duplicated reads (GATK best practice workflow)
#---
#This script automates the marking of duplicate reads using GATKs MarkDuplicates to create
#analysis-ready reads for variant call analysis.

#define global variables
DIR="/path/to/directory" #working directory

#create required directories
mkdir -p "${DIR}/outputs/analysis-ready_reads"


#Marking duplicate reads (markDuplicates)
#---
#create an array of sample names from folder names
mapfile -t sample_names < <(find ${DIR}/samples -maxdepth 1 -type d ! -path ${DIR}/samples -exec basename {} \;)

#create empty arrays
valid_files=()
invalid_files=()
for sample_name in ${sample_names[@]}; do
	#record start time
	start_time=$(date +%s)

	#define the input sorted BAM file
	input_bam="${DIR}/outputs/alignments/${sample_name}.bam"

	#define the output BAM file with analysis-ready reads
	output_bam="${DIR}/outputs/analysis-ready_reads/${sample_name}.bam"

	#define the metrics file
	metrics_file="${DIR}/outputs/analysis-ready_reads/dup_metrics-${sample_name}.txt"

	#run the duplicate flagger
	echo "Starting duplicate read marking for sample ${sample_name}..."
	gatk MarkDuplicates -I $input_bam -O $output_bam -M $metrics_file \
	2> "${DIR}/logs/alignments/MarkDuplicates-${sample_name}.log"

	#logic to deal with failed files
	if [ $? -eq 0 ]; then
		#update valid list of files
		echo "Sample ${sample_name} was successful!"
		valid_files+=("${sample_name}")
	else
		#update the invalid list of files
		echo "Sample ${sample_name} failed!"
		invalid_files+=("${sample_name}")
	fi

	#record the end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time))
	minutes=$(( run_time / 60 ))
    echo "Sample ${sample_name} took ${minutes} minutes"
done

#print the failed files
echo "Duplicate read marking of BAM files has completed!"
echo "The following files failed: ${invalid_files[@]}"
