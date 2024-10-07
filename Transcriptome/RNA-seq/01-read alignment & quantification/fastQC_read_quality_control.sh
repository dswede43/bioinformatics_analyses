#!/bin/bash

#Read quality control - FastQC
#---
#This script automates the quality control of all FASTQ files using FastQC.

#define global variables
DIR="/path/to/directory" #working directory

#create directories
mkdir -p "${DIR}/outputs/quality control"
mkdir -p "${DIR}/logs/quality control"


#Sample quality control
#---
#create an array of sample names from folder names
mapfile -t sample_names < <(find ${DIR}/samples -maxdepth 1 -type d ! -path ${DIR}/samples -exec basename {} \;)

#create empty arrays
run_times=()
valid_files=()
invalid_files=()
for sample_name in ${sample_names[@]}; do
	#record start time
	start_time=$(date +%s)

	#input samples
	mapfile -t inputs \
    < <(find samples/${sample_name} -type f -name "*.fastq.gz")

	input_num=1
	for input in ${inputs[@]}; do
		echo "Starting quality control for sample ${sample_name}..."

		#run the quality control
		fastqc "$DIR/${input}" -o "$DIR/outputs/quality control" \
		2> "${DIR}/logs/quality control/fastqc-${sample_name}_${input_num}.log"

		#logic to deal with failed files
		if [ $? -eq 0 ]; then
			#update valid list of files
			echo "Sample ${sample_name} was successful!"
			valid_files+=("${sample_name}-${input_num}")
		else
			#update the invalid list of files
			echo "Sample ${sample_name} failed!"
			invalid_files+=("${sample_name}-${input_num}")
		fi

		#increase the input number
		input_num=$((input_num+1))
	done

	#record the end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time))
	minutes=$(( run_time / 60 ))
	run_times+=("sample ${sample_name}: ${minutes} minutes")
    echo "Sample ${sample_name} took ${minutes} minutes"
done

#print the failed SRA files
echo "Quality control of FASTQ files has completed!"
echo "The following files failed: ${invalid_files[@]}"
