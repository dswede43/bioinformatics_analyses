#!/bin/bash

#Read quantification-featureCounts
#---
#This script automates the quantification of aligned reads from SAM files using featureCounts.
#featureCounts documentation: https://subread.sourceforge.net/SubreadUsersGuide.pdf

#directory structure
#---samples/
#---outputs/quants/
#---annotations/

#define global variables
DIR="/media/sf_VM_share/RNA-seq_data/KneeOA" #working directory


#Read quantification
#---
#define the genome annotation file
genome_annotation="$DIR/annotation/Homo_sapiens.GRCh38.110.gtf"

#create an array of sample names from folder names
mapfile -t sample_names < <(find "$DIR/samples" -maxdepth 1 -type d ! -path "$DIR/samples" -exec basename {} \;)

#create an empty list of run times
run_times=()

#for each sample
for sample_name in ${sample_names[@]}; do

    #record start time
    start_time=$(date +%s)

    #alignment file
    mapfile -t alignment_file \
    < <(find "$DIR/outputs/mappings/$sample_name" -type f -name "*.sam")

    #create an output file
    mkdir -p $DIR/outputs/quants/$sample_name
    output_file="$DIR/outputs/quants/$sample_name/sample$sample_name.txt"

    #run the quantification
    echo "Starting quantification of sample $sample_name..."
    featureCounts -s 2 -T 8 -a $genome_annotation -o $output_file $alignment_file 2> "$DIR/outputs/quants/$sample_name/quantification.log"
    echo "Quantification of sample $sample_name has finished!"

	#record the end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time ))
	run_times+=("sample $sample_name: $(( run_time / 60 )) minutes")
    echo "The quantification for sample $sample_name took "$(( run_time / 60 ))" minutes"
done

#save the sample run times as .txt file
printf "%s\n" "${run_times[@]}" > "$DIR/outputs/quants/quantification_run_times.txt"
