#!/bin/bash

#Read alignmen - BWA (Burrow-Wheeler Aligner)
#---
#This script automates the alignment of reads to a reference genome using BWA-MEM.


#Table of contents:
#---
#1. Reference genome indexing
#2. Read alignment to reference genome
#3. Compress outputs to .BAM and sort by coordinate

#define global variables
DIR="/path/to/directory" #working directory
PAIRED=true #paired-end reads?

#create required directories
mkdir -p "$DIR/outputs/alignments"
mkdir -p "$DIR/logs/read alignment"


#1. Reference genome indexing
#---
#define the reference index FASTA file
ref="${DIR}/OmykA1.1/GCF_013265735.2/GCF_013265735.2_USDA_OmykA_1.1_genomic.fna"

if false; then
	#build the reference index if it has not been built yet
	echo "Building the reference index..."
	bwa index ${ref}
	echo "Reference index complete!"
fi


#2. Read alignment to reference genome
#---
#create an array of sample names from folder names
mapfile -t sample_names < <(find "${DIR}/samples" -maxdepth 1 -type d ! -path "${DIR}/samples" -exec basename {} \;)

#create empty arrays
run_times=()
valid_files=()
invalid_files=()
for sample_name in ${sample_names[@]}; do
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
		echo "Starting alignment for sample ${sample_name}..."
		bwa mem -t 8 -R "@RG\tID:${sample_name}\tPL:ILLUMINA\tSM:${sample_name}" $ref $input_1 $input_2 > $output \
		2> "${DIR}/logs/read alignment/bwa-mem-${sample_name}.log"

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
		echo "Starting alignment for sample ${sample_name}..."
		bwa mem -t 8 -R "@RG\tID:${sample_name}\tPL:ILLUMINA\tSM:${sample_name}" $ref $input > $output \
		2> "${DIR}/logs/read alignment/bwa-mem-${sample_name}.log"

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
	fi

	#record end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time))
	minutes=$(( run_time / 60 ))
    echo "Sample ${sample_name} took ${minutes} minutes"
done

#print the failed samples
echo "Read alignment has completed!"
echo "The following samples failed: ${invalid_files[@]}"


#3. Compress outputs to .BAM and sort by coordinate
#---
for sample_name in ${sample_names[@]}; do
	#define the input SAM file
	sam_file="${DIR}/outputs/alignments/${sample_name}.sam"

	#define the output BAM file
	bam_file="${DIR}/outputs/alignments/${sample_name}.bam"

	#convert .SAM to .BAM
	echo "Compressing to ${sample_name} to .BAM..."
	samtools view -Sb $sam_file > $bam_file
	echo "Compression to .BAM complete!"

	#sort the .BAM file
	sorted_bam_file="${DIR}/outputs/alignments/${sample_name}_sorted.bam"
	echo "Sorting coordinates for ${sample_name}..."
	samtools sort $bam_file -o $sorted_bam_file
	echo "Sorting complete!"
done
