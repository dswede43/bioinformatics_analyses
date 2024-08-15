#!/bin/bash

#Part 1: TIGeR bioinformatics interview assignment
#---
#This script automates various workflows in GATK to filter VCF files.


#Table of contents:
#---
#1. Accept multiple vcf files as inputs
#2. Perform integrity checks
#3. Hard-filter variants

#define global variables
DIR="/gatk/my_data"


#1. Accept multiple vcf files as inputs
#---
#define an empty array
vcf_files=()

for vcf_file in "$DIR"/*.vcf.gz; do
    if [ -f "$vcf_file" ]; then
        #append the vcf file name to the array
        vcf_files+=("$vcf_file")
    fi
done


#2. Perform integrity checks
#---
echo "Validating integrity of VCF files..."

#define empty arrays
valid_vcf_files=()
invalid_vcf_files=()
for vcf_file in "${vcf_files[@]}"; do
    #create the index file for the current VCF file
    gatk IndexFeatureFile \
        -F "$vcf_file" \
        2> "$DIR/logs/IndexFeatureFile_logs.log"

    #validate the integrity of the current VCF file
    gatk ValidateVariants \
        -V "$vcf_file" \
        2> "$DIR/logs/ValidateVariants_logs.log"

    #logic to deal with invalid VCF files
    if [ $? -eq 0 ]; then
        #update the valid list of VCF files
        echo "${vcf_file##*/} is valid."
        valid_vcf_files+=("$vcf_file")

    else
        #update the invalid list of VCF files
        echo "${vcf_file##*/} is invalid."
        invalid_vcf_files+=("$vcf_file")
    fi
done

#update the original input of VCF files
vcf_files=("${valid_vcf_files[@]}")

#print the list of invalid VCF files
echo "Invalid VCF files that have been removed:"
for vcf_file in "${invalid_vcf_files[@]}"; do
    echo "${vcf_file##*/}"
done


#3. Hard-filter variants
#---
echo "Applying hard-filtering to variants for remaining VCF files..."

for vcf_file in "${vcf_files[@]}"; do
    #define the output VCF file path
    vcf_output="$DIR/filtered_vcfs/${vcf_file##*/}"

    #subset to SNP variants only
    gatk SelectVariants \
        -V "$vcf_file" \
        -select-type SNP \
        -O "${vcf_output%.vcf.gz}_snps.vcf.gz" \
        2> "$DIR/logs/SelectVariants_logs.log"

    #apply hard-filtering to variants
    gatk VariantFiltration \
        -V "${vcf_output%.vcf.gz}_snps.vcf.gz" \
        -filter "QUAL < 20" --filter-name "QUAL20" \
        -filter "DP < 5" --filter-name "DP5" \
        -O "${vcf_output%.vcf.gz}_snps_filtered.vcf.gz" \
        2> "$DIR/logs/VariantFiltration_logs.log"

    #remove input VCF file
    rm -rf "${vcf_output%.vcf.gz}_snps.vcf.gz"
    rm -rf "${vcf_output%.vcf.gz}_snps.vcf.gz.tbi"
done

echo "Variant hard-filtering complete! VCF output files found in $DIR/filtered_vcfs."

