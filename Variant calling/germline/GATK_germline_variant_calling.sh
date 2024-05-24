#!/bin/bash

#Variant calling for germline mutations
#---
#This script identifies germline basepair variations including single nucleotide polymorphisms (SNPs),
#indels, and structural variations (duplications, inversions, translocations) via HalpotypeCaller.

#dependencies:
#broadinstitute/gatk docker image (https://hub.docker.com/r/broadinstitute/gatk)

#directory structure:
#---references
#---outputs/mappings
#---outputs/variants

#define global variables
DIR="/media/sf_VM_share/RNA-seq_data/1000_genomes_project"

#start the GATK docker container
docker run -v ${DIR}:/gatk/my_data -it broadinstitute/gatk:latest

#change directory
cd my_data


#Germline variant caller (HaplotypecCaller)
#---
#create an array of sample names from folder names
mapfile -t sample_names < <(find samples -maxdepth 1 -type d ! -path samples -exec basename {} \;)

#define reference index file
ref_index="references/hg38.fa"

#create an empty list of run times
run_times=()

for sample_name in ${sample_names}; do
	#record start time
	start_time=$(date +%s)

	#define input BAM file with analysis-ready reads
	input_bam=outputs/mappings/${sample_name}/sample${sample_name}_dup_flags_bqsr.bam

	#define the output VCF file
    mkdir -p outputs/variants/$sample_name
	output_vcf=outputs/mappings/${sample_name}/sample${sample_name}_raw_variants.vcf

	#run germline variant calling
	gatk HaplotypeCaller -I $input_bam -O $output_vcf -R $ref_index

	#record the end time
	end_time=$(date +%s)

	#store the sample runtime
	run_time=$(( end_time - start_time ))
	run_times+=("sample $sample_name: $(( run_time / 60 )) minutes")
    echo "Base pair quality recalibration for sample $sample_name took "$(( run_time / 60 ))" minutes"
done

#save the sample run times as .txt file
printf "%s\n" "${run_times[@]}" > "$DIR/outputs/mappings/runtimes/variant_call_run_times.txt"
