#!/bin/bash

#Read quality control-FastQC
#---
#This script automates the quality control of all FASTQ files using FastQC.

#directory structure:
#---samples/

#define global variables
DIR="/path/to/directory" #working directory


#Sample quality control
#---
#create an array of sample names from folder names
mapfile -t sample_names < <(find ${DIR}/samples -maxdepth 1 -type d ! -path ${DIR}/samples -exec basename {} \;)

#create an empty list of run times
run_times=()

for sample_name in ${sample_names}; do
	#record start time
	start_time=$(date +%s)

	#input samples
	mapfile -t inputs \
    < <(find samples/${sample_name} -type f -name "*.fastq.gz")

	for input in ${inputs}; do
		#run the quality control
		fastqc ${input} -o samples/${sample_name}
	done

	#record the end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time ))
	run_times+=("sample $sample_name: $(( run_time / 60 )) minutes")
    echo "The read quality control for sample $sample_name took "$(( run_time / 60 ))" minutes"
done
