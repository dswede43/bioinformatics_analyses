#!/bin/bash

#Paired-end read alignment with HISAT2
#---
#This script automates the alignment of paired-end reads from FASTQ files to a
#reference index using HISAT2.


#define global variables
DIR="/path/to/directory" #working directory
REF_INDEX="$DIR/grch38/genome" #reference index directory

#define the range of SRA sample id's
SRA_MIN=6298258
SRA_MAX=6298287

#create required directories
mkdir -p "$DIR/outputs/alignments"
mkdir -p "$DIR/logs/alignments"


#Paired-end read alignment
#---
#create a list of SRA sample id's
echo "Creating SRA identity names..."

sample_names=()
for ((i = SRA_MIN; i <= SRA_MAX; i++)); do
    sample_names+=("SRR$i")
done

#create empty arrays
run_times=()
valid_alignments=()
invalid_alignments=()
for sample_name in ${sample_names[@]}; do
    #record start time
    start_time=$(date +%s)

    #prepare the fastq paried-end input lists
    mapfile -t input_list \
    < <(find "$DIR/samples/fastq" -type f -name "*$sample_name*")

    input_list1=""
    input_list2=""
    for fastq_file in ${input_list[@]}; do
        if [[ "$fastq_file" == *"_1"* ]]; then
            input_list1+=",$fastq_file"
        else
            input_list2+=",$fastq_file"
        fi
    done

    #remove the initial commas
    input_list1="${input_list1#,}"
    input_list2="${input_list2#,}"

    #create an output file
    output="$DIR/outputs/alignments/$sample_name.sam"

    #run the alignment
    echo "Starting alignment of sample $sample_name..."
    hisat2 -p 8 -q --rna-strandness R -x $REF_INDEX -1 $input_list1 -2 $input_list2 -S $output 2> "$DIR/logs/alignments/${sample_name}_alignment.log"

    #logic to deal with failed alignments
    if [ $? -eq 0 ]; then
        #update valid list of aligned samples
        echo "${sample_name} alignment completed successfully!"
        valid_alignments+=("$sample_name")

        #delete fastq files to save space
        find "$DIR/samples/fastq" -type f -name "*$sample_name*" -exec rm -f {} \;
    else
        #update invalid list of aligned samples
        echo "${sample_name} alignment failed!"
        invalid_alignments+=("$sample_name")
    fi

	#record the end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time ))
	run_times+=("$sample_name alignment runtime: $(( run_time / 60 )) minutes")
    echo "$sample_name alignment runtime: "$(( run_time / 60 ))" minutes"
done

#print the failed alignments
echo "Read alignments complete!"
echo "The following FASTQ files failed to align: ${invalid_alignments[@]}"

#save the sample run times as .txt file
printf "%s\n" "${run_times[@]}" > "$DIR/outputs/alignment_run_times.txt"