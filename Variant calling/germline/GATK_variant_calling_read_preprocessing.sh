#!/bin/bash

#Read preprocessing for variant calling analysis of germline mutations (GATK best practice workflow)
#---

#This script completes the preprocessing steps from raw unmapped reads (FASTQ) to
#analysis-ready reads (BAM) required for variant calling analysis of germline mutations
#using GATK best practice workflow. Steps include read duplicate flagging (markDuplicateSpark)
#and base pair quality score recalibration.

#dependencies:
#broadinstitute/gatk docker image (https://hub.docker.com/r/broadinstitute/gatk)

#directory structure:
#---references
#---outputs/mappings

#define global variables
DIR="/media/sf_VM_share/RNA-seq_data/1000_genomes_project"

#start the GATK docker container
docker run -v ${DIR}:/gatk/my_data -it broadinstitute/gatk:latest

#change directory
cd my_data


#Flag duplicate reads
#---
#create an array of sample names from folder names
mapfile -t sample_names < <(find samples -maxdepth 1 -type d ! -path samples -exec basename {} \;)

#create an empty list of run times
dup_flag_run_times=()

for sample_name in ${sample_names}; do
	#record start time
	start_time=$(date +%s)


	#Read sorting
	#---
	#define the input SAM file
	input_sam=outputs/mappings/${sample_name}/sample${sample_name}.sam

	#define the output BAM file
	output_bam=outputs/mappings/${sample_name}/sample${sample_name}_sorted.bam

	#sort reads by alignment positions (coordinates)
	gatk SortSam -I $input_sam -O $output_bam -SO coordinate


	#Read duplicate flagging
	#---
	#define the input BAM file
	input_bam=outputs/mappings/${sample_name}/sample${sample_name}_sorted.bam

	#define the output BAM file
	output_bam=outputs/mappings/${sample_name}/sample${sample_name}_dup_flags.bam

	#run the duplicate flagger
	gatk MarkDuplicates -I $input_bam -O $output_bam -M outputs/mappings/${sample_name}/dup_metrics.txt

	#record the end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time ))
	dup_flag_run_times+=("sample $sample_name: $(( run_time / 60 )) minutes")
    echo "The duplicate read flagger for sample $sample_name took "$(( run_time / 60 ))" minutes"
done

#save the sample run times as .txt file
printf "%s\n" "${dup_flag_run_times[@]}" > "$DIR/outputs/mappings/runtimes/dup_flag_run_times.txt"


#Base pair quality score recalibration
#---
#define the VCF with known variant sites
known_sites="references/Homo_sapiens_assembly38.dbsnp138.vcf"

#define reference index file
ref_index="references/hg38.fa"

#create an empty list of run times
bqsr_run_times=()

for sample_name in ${sample_names}; do
	#record start time
	start_time=$(date +%s)


	#Recalibration table
	#---
	#define the input BAM file with flagged duplicate reads
	input_bam=outputs/mappings/${sample_name}/sample${sample_name}_dup_flags.bam

	#define the output recalibration table file
	recal_table=outputs/mappings/${sample_name}/sample${sample_name}_recal_data.table

	#generate the recalibration table based on various covariates
	gatk BaseRecalibrator -I $input_bam --known-sites $known_sites -O $recal_table -R $ref_index


	#Apply recalibration
	#---
	#define the output BAM file with recalibrated base scores
	output_bqsr=outputs/mappings/${sample_name}/sample${sample_name}_dup_flags_bqsr.bam

	#apply the model to adjust the base quality scores
	gatk ApplyBQSR --bqsr-recal-file $recal_table -I $input_bam -O $output_bqsr

	#record the end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time ))
	bqsr_flag_run_times+=("sample $sample_name: $(( run_time / 60 )) minutes")
    echo "Base pair quality recalibration for sample $sample_name took "$(( run_time / 60 ))" minutes"
	
	#remove unnecessary files to save space
	files=("outputs/mappings/${sample_name}/*.sam" "outputs/mappings/${sample_name}/*sorted.bam" "outputs/mappings/${sample_name}/*dup_flags.bam")
	for file in ${files[@]}; do
		rm -v $file
	done
done

#save the sample run times as .txt file
printf "%s\n" "${bqsr_run_times[@]}" > "$DIR/outputs/mappings/runtimes/bqsr_run_times.txt"
