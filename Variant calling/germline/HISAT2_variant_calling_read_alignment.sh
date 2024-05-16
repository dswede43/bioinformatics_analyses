#!/bin/bash

#Read alignment for variant calling analysis of germline mutations (GATK best practice workflow)
#---
#This script completes the paired-end read alignment using HISAT2 to create analysis-ready reads (BAM)
#required for variant calling analysis of germline mutations using GATK best practice workflow.

#dependencies:
#broadinstitute/gatk docker image (https://hub.docker.com/r/broadinstitute/gatk)

#directory structure:
#---samples
#---references
#---outputs/mappings

#define global variables
DIR="/media/sf_VM_share/RNA-seq_data/1000_genomes_project"


#Reference indexing preparation for HISAT2
#---
if false; then
	#build index from FASTA file for HISAT2
	hisat2-build -p 8 "$DIR/references/gr38.fa" genome
fi


#Read alignment using HISAT2
#---
#create an empty list of run times
run_times=()

#define the reference genome index directory
ref_index="$DIR/references/grch38/genome"

#create an array of sample names from folder names
mapfile -t sample_names < <(find samples -maxdepth 1 -type d ! -path samples -exec basename {} \;)

for sample_name in ${sample_names}; do
	#record start time
	start_time=$(date +%s)

	#input lists of paired FASTQ files
    #mate 1
    mapfile -t input_list1 \
    < <(find "$DIR/samples/$sample_name" -type f -name "*1.fastq*")
    #mate 2
    mapfile -t input_list2 \
    < <(find "$DIR/samples/$sample_name" -type f -name "*2.fastq*")

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
    hisat2 -p 8 -q -x $ref_index -1 $input_1 -2 $input_2 -S $output 2> "$DIR/outputs/mappings/$sample_name/alignment.log"
    echo "Alignment of sample $sample_name finished!"

	#record the end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time ))
	run_times+=("sample $sample_name: $(( run_time / 60 )) minutes")
    echo "The alignment for sample $sample_name took "$(( run_time / 60 ))" minutes"
done

#save the sample run times as .txt file
printf "%s\n" "${run_times[@]}" > "$DIR/outputs/mappings/runtimes/alignment_run_times.txt"


#Adding read group information
#---
#start the GATK docker container
docker run -v ${DIR}:/gatk/my_data -it broadinstitute/gatk:latest

#change directory
cd my_data

#create an array of sample names from folder names
mapfile -t sample_names < <(find samples -maxdepth 1 -type d ! -path samples -exec basename {} \;)

for sample_name in ${sample_names}; do
	#define the input SAM file
	input_sam=outputs/mappings/${sample_name}/sample${sample_name}.sam

	#define the output SAM file
	output_sam=outputs/mappings/${sample_name}/sample${sample_name}_RG.sam

	#add read group information
	gatk AddOrReplaceReadGroups -I $input_sam -O $output_sam \
		-ID $sample_name \
		-LB lib1 \
		-PL ILLUMINA \
		-PU unit1 \
		-SM $sample_name
done
