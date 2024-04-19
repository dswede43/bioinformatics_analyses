#!/bin/bash

#Paired-end read alignment-HISAT2
#---
#This script automates the alignment of paired-end reads from FASTQ files for multiple samples
#to a reference index and output SAM files.

#directory structure:
#---samples/
#---<reference_index_name>/
#---outputs/mappings/

#define global variables
DIR="/path/to/directory" #working directory


#Read mapping
#---
#define the genome index directory
reference_index="$DIR/grch38/genome"

#create an array of sample names from folder names
mapfile -t sample_names < <(find "$DIR/samples" -maxdepth 1 -type d ! -path "$DIR/samples" -exec basename {} \;)

#create an empty list of run times
run_times=()

#for each sample
for sample_name in ${sample_names[@]}; do

	#record start time
	start_time=$(date +%s)

    #input lists of paired FASTQ files
    #mate 1
    mapfile -t input_list1 \
    < <(find "$DIR/samples/$sample_name" -type f -name "*1.fastq")
    #mate 2
    mapfile -t input_list2 \
    < <(find "$DIR/samples/$sample_name" -type f -name "*2.fastq")

    #convert the input list to a single comma-separated string
    #mate 1
    input_1=""
    for i in "${input_list1[@]}"; do
        input_1+=",$i"
    done
    input_1="${input_1#,}" #remove the initial comma from the string
    #mate 2
    input_2=""
    for i in "${input_list2[@]}"; do
        input_2+=",$i"
    done
    input_2="${input_2#,}" #remove the initial comma from the string

    #create an output file
    mkdir -p $DIR/outputs/mappings/$sample_name
    output="$DIR/outputs/mappings/$sample_name/sample$sample_name.sam"

	#run the alignment
    echo "Starting alignment of sample $sample_name..."
    hisat2 -p 8 -q --rna-strandness R -x $reference_index -1 $input_1 -2 $input_2 -S $output 2> "$DIR/outputs/mappings/$sample_name/alignment.log"
    echo "Alignment of sample $sample_name finished!"

	#record the end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time ))
	run_times+=("sample $sample_name: $(( run_time / 60 )) minutes")
    echo "The alignment for sample $sample_name took "$(( run_time / 60 ))" minutes"
done

#save the sample run times as .txt file
printf "%s\n" "${run_times[@]}" > "$DIR/outputs/mappings/alignment_run_times.txt"
