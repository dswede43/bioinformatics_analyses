#!/bin/bash

#Read alignment with HISAT2
#---
#This script automates the alignment of reads from FASTQ files to a
#reference index using HISAT2.


#define global variables
DIR="/path/to/directory" #working directory
REF_INDEX="$DIR/grch38/genome" #reference index directory

#create required directories
mkdir -p "$DIR/outputs/alignments"
mkdir -p "$DIR/logs/alignments"


#Read alignment
#---
#create an array of sample names from folder names
mapfile -t sample_names < <(find samples -maxdepth 1 -type d ! -path samples -exec basename {} \;)

#create empty arrays
valid_files=()
invalid_files=()
for sample_name in ${sample_names}; do
	#record start time
	start_time=$(date +%s)

	#paired-end input samples
	if [ $PAIRED == true ]; then
		#mate 1
		mapfile -t input_list1 \
		< <(find "${DIR}/samples/${sample_name}" -type f -name "*1.fastq*")

		#convert the input list into single comma-separated string
		input_1=""
		for i in "${input_list1[@]}"; do
			input_1+=",$i"
		done
		input_1="${input_1#,}" #remove the initial comma from the string

		#mate 2
		mapfile -t input_list2 \
		< <(find "${DIR}/samples/${sample_name}" -type f -name "*2.fastq*")

		#convert the input list into single comma-separated string
		input_2=""
		for i in "${input_list2[@]}"; do
			input_2+=",$i"
		done
		input_2="${input_2#,}" #remove the initial comma from the string

		#define output file
		output="${DIR}/outputs/alignments/${sample_name}.sam"

		#run the alignment
		echo "Starting alignment of sample $sample_name..."
		hisat2 -p 8 -q --rna-strandness R -x $ref_index -1 $input_1 -2 $input_2 -S $output \
		2> "${DIR}/logs/read alignment/HISAT2-${sample_name}.log"

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

	#single-end input samples
	else
		mapfile -t input_list \
		< <(find "${DIR}/samples/${sample_name}" -type f -name "*.fastq*")

		#convert the input list into single comma-separated string
		input=""
		for i in "${input_list[@]}"; do
			input+=",$i"
		done
		input="${input#,}" #remove the initial comma from the string

		#define output file
		output="${DIR}/outputs/alignments/${sample_name}.sam"

		#run the alignment
		echo "Starting alignment of sample $sample_name..."
		hisat2 -p 8 -q --rna-strandness R -x $REF_INDEX -U $input -S $output \
		2> "${DIR}/logs/read alignment/HISAT2-${sample_name}.log"
	fi

	#record the end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time))
	minutes=$(( run_time / 60 ))
    echo "Sample ${sample_name} took ${minutes} minutes"
done

#print the failed samples
echo "Read alignment has completed!"
echo "The following samples failed: ${invalid_files[@]}"
