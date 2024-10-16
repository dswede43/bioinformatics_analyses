#!/bin/bash

#Read base pair quality score recalibration (BQSR)
#---
#Script to automate the recalibration of read base pair quality scores as part of
#GATKs best practice workflow.

#define global variables
DIR="/path/to/directory"

#create required directories
mkdir -p "$DIR/outputs/analysis-ready_reads"
mkdir -p "$DIR/logs/alignments"


#Base pair quality score recalibration (BQSR)
#---
#define the VCF with known variant sites
known_sites="references/known_sites.vcf"

#define reference index file
ref_index="references/ref_index.fa"

#create an array of sample names from folder names
mapfile -t sample_names < <(find ${DIR}/samples -maxdepth 1 -type d ! -path ${DIR}/samples -exec basename {} \;)

#create empty arrays
valid_files=()
invalid_files=()
for sample_name in ${sample_names[@]}; do
	#record start time
	start_time=$(date +%s)


	#1. Recalibration table (BaseRecalibrator)
	#---
	#define the input BAM file with flagged duplicate reads
	input_bam="${DIR}/outputs/analysis-ready_reads/${sample_name}.bam"

	#define the output recalibration table file
	recal_table="${DIR}/outputs/analysis-ready_reads/recal_table-${sample_name}.table"

	#generate the recalibration table based on various covariates
	echo "Generating the recalibration table for sample ${sample_name}..."
	gatk BaseRecalibrator -I $input_bam --known-sites $known_sites -O $recal_table -R $ref_index \
	2> "${DIR}/logs/alignments/BaseRecalibrator-${sample_name}.log"


	#2. Apply recalibration (ApplyBQSR)
	#---
	#define the output BAM file with recalibrated base scores
	output_bqsr="${DIR}/outputs/analysis-read_reads/${sample_name}.bam"

	#apply the model to adjust the base quality scores
	echo "Adjust the base quality scores for sample ${sample_name}..."
	gatk ApplyBQSR --bqsr-recal-file $recal_table -I $input_bam -O $output_bqsr \
	2> "${DIR}/logs/alignments/ApplyBQSR-${sample_name}.log"
	echo "Base pair recalibration has completed!"

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
echo "Base pair recalibration has completed!"
echo "The following files failed: ${invalid_files[@]}"
