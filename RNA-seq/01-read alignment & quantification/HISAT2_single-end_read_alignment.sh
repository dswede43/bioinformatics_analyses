#!/bin/bash

#Single-end read alignment-HISAT2
#---
#This script automates the alignment of single-end reads from FASTQ files for multiple samples
#to a reference index and output SAM files.
#HSIAT2 documentation: https://daehwankimlab.github.io/hisat2/manual/

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

    #input list of FASTQ files
    mapfile -t input_list \
    < <(find "$DIR/samples/$sample_name" -type f -name "*.fastq.gz")

    #convert the input list to a single comma-separated string
    input=""
    for i in "${input_list[@]}"; do
        input+=",$i"
    done
    input="${input#,}" #remove the initial comma from the string

    #create an output file
    mkdir -p $DIR/outputs/mappings/$sample_name
    output="$DIR/outputs/mappings/$sample_name/sample$sample_name.sam"

	#run the alignment
    echo "Starting alignment of sample $sample_name..."
    hisat2 -p 8 -q --rna-strandness R -x $reference_index -U $input -S $output 2> "$DIR/outputs/mappings/$sample_name/alignment.log"
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
