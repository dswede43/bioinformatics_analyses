#!/bin/bash

#Batch downloading paired-end FASTQ files
#---
#This script automates the batch downloading of SRA files and extraction
#into FASTQ files containing paired-end reads using the SRA toolkit.
#SRA-toolkits documentation: https://github.com/ncbi/sra-tools/wiki


#Table of contents:
#---
#1. Download SRA files
#2. Convert SRA files into FASTQ files
#3. Compress all FASTQ files using gzip

#define global variables
DIR="/path/to/directory" #working directory

#define the range of SRA sample id's
SRA_MIN=6298258
SRA_MAX=6298307

#create required directories
mkdir -p "$DIR/samples/fastq"
mkdir -p "$DIR/samples/sra"


#1. Download SRA files
#---
#create a list of SRA sample id's
echo "Creating SRA identity names..."

sra_ids=()
for ((i = SRA_MIN; i <= SRA_MAX; i++)); do
    sra_ids+=("SRR$i")
done

#download sample SRA files
echo "Downloading SRA files..."

#define empty arrays
valid_sra_ids=()
invalid_sra_ids=()
for sra_id in ${sra_ids[@]}; do
    #download the current SRA file
    echo "Downloading SRA file for $sra_id..."
    prefetch "$sra_id"
        --progress \ #display progress bar
        --output-directory "$DIR/samples" \ #output directory
        2> "$DIR/logs/prefetch_logs.log"

    #logic to deal with failed downloads
    if [ $? -eq 0 ]; then
        #update valid lists of SRA files
        echo "${sra_id} was successfully downloaded!"
        valid_sra_ids+=("$sra_id")
    else
        #update the invalid list of SRA files
        echo "${sra_id} failed to download!"
        invalid_sra_ids+=("$sra_id")
    fi
done

#print the failed SRA files
echo "SRA file downloads complete!"
echo "The following SRA files failed to download: ${invalid_sra_ids[@]}"


#2. Convert SRA files into FASTQ files
#---
#define empty arrays
valid_sra_dumps=()
invalid_sra_dumps=()
for sra_id in ${sra_ids[@]}; do
    #convert the current SRA file to FASTQ
    echo "Converting to FASTQ for $sra_id..."

    #split-3: splits paired-end reads into seperate fastq files
    #skip-technical: ignores technical reads (adapters, primers, barcodes, etc)
    fasterq-dump "$sra_id" --split-3 --skip-technical --progress --outdir "$DIR/samples/fastq" \
        2> "$DIR/logs/fasterq-dump_logs.log"

    #logic to deal with failed conversions
    if [ $? -eq 0 ]; then
        #update valid lists of SRA dumps
        echo "${sra_id} was successfully converted into FASTQ!"
        valid_sra_dumps+=("$sra_id")
    else
        #update the invalid list of SRA dumps
        echo "${sra_id} failed to convert into FASTQ!"
        invalid_sra_dumps+=("$sra_id")
    fi

    #move the SRA file
    mv "$DIR/samples/$sra_id/"*.sra "$DIR/samples/sra/"
    rm -rf "$DIR/samples/$sra_id"
done

#print the failed FASTQ conversions
echo "FASTQ conversions complete!"
echo "The following SRA files failed to convert to FASTQ: ${invalid_sra_dumps[@]}"


#3. Compress all FASTQ files using gzip
gzip "$DIR/samples/fastq/"*.fastq
