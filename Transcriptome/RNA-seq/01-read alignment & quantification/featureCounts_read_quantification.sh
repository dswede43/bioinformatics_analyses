#!/bin/bash

#Read quantification with featureCounts
#---
#This script automates the quantification of aligned reads from SAM files
#using featureCounts.
#featureCounts documentation: https://subread.sourceforge.net/SubreadUsersGuide.pdf


#define global variables
DIR="/media/sf_VM_share/RNA-seq_data/Immune_cells-GSE107011" #working directory
GENOME_ANNOTATION="$DIR/annotation/Homo_sapiens.GRCh38.110.gtf" #genome annotation file

#define the range of SRA sample id's
SRA_MIN=6298258
SRA_MAX=6298287

#create required directories
mkdir -p "$DIR/outputs/quants"
mkdir -p "$DIR/logs/quants"


#Read quantification
#---
#create a list of SRA sample id's
echo "Creating SRA identity names..."

sample_names=()
for ((i = SRA_MIN; i <= SRA_MAX; i++)); do
    sample_names+=("SRR$i")
done

#create empty arrays
run_times=()
valid_quants=()
invalid_quants=()
for sample_name in ${sample_names[@]}; do
    #record start time
    start_time=$(date +%s)

    #prepare the input SAM file
    mapfile -t sam_file \
    < <(find "$DIR/outputs/alignments" -type f -name "*$sample_name*")

    #create an output file
    output="$DIR/outputs/quants/$sample_name.txt"

    #run the quantification
    echo "Starting quantification of sample $sample_name..."
    featureCounts -s 2 -T 8 -a $GENOME_ANNOTATION -o $output_file $sam_file \
    2> "$DIR/logs/quants.log"

    #logic to deal with failed quantifications
    if [ $? -eq 0 ]; then
        #update valid list of quantified samples
        echo "${sample_name} quantification completed successfully!"
        valid_quants+=("$sample_name")

    else
        #update invalid list of quantified samples
        echo "${sample_name} quantification failed!"
        invalid_quants+=("$sample_name")
    fi

	#record the end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time ))
	run_times+=("$sample_name alignment runtime: $(( run_time / 60 )) minutes")
    echo "$sample_name alignment runtime: "$(( run_time / 60 ))" minutes"
done

#print the failed quantifications
echo "Read quantification complete!"
echo "The following SAM files failed to be quantifies: ${invalid_quants[@]}"

#save the sample run times as .txt file
printf "%s\n" "${run_times[@]}" > "$DIR/outputs/quant_run_times.txt"
