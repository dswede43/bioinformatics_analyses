#!/bin/bash

#Batch downloading FASTQ files
#---
#This script automates the batch downloading of SRA files and extraction
#into FASTQ files containing paired-end reads using the SRA toolkit.
#SRA-toolkits documentation: https://github.com/ncbi/sra-tools/wiki

#required directory structure:
#---samples/

#define global variables
DIR="/path/to/directory" #working directory

#define the range of SRA sample id's
SRA_MIN=6298258
SRA_MAX=6298259


#Download FASTQ files
#---
#create a list of SRA sample id's
sra_ids=()
for ((i = SRA_MIN; i <= SRA_MAX; i++)); do
    sra_ids+=("SRR$i")
done

#download sample SRA files
for sra_id in ${sra_ids[@]}; do
    echo "Downloading SRA file for $sra_id..."
    prefetch "$sra_id" --progress --output-directory "$DIR/samples"
done

#Convert SRA files into FASTQ files
#---
#--split-3: splits paired-end reads into *_1.fastq and *_2.fastq files, unmated reads are placed into *.fastq
#--skip-technical: ignores technical reads (adapters, primers, barcodes, etc)
for sra_id in ${sra_ids[@]}; do
    echo "Extracting FASTQ file for $SRA_ID..."
    fasterq-dump "$SRA_ID" --progress --split-3 --skip-technical --outdir "$DIR/samples/$sra_id"

    #remove the .sra file
    echo "Deleting SRA file to save space"
    rm -rfv "$DIR/samples/$SRA_ID/*.sra"
done
