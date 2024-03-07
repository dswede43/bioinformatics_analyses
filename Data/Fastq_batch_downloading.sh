#!/bin/bash

#Batch downloading FASTQ files
#---
#This script automates batch downloading of FASTQ files using the SRA toolkit
#(https://github.com/ncbi/sra-tools/wiki)

#define SRA sample id's
SRA_MIN=6298258
SRA_MAX=6298384

SRA_IDS=()
for ((i = SRA_MIN; i <= SRA_MAX; i++)); do
    SRA_IDS+=("SSR$i")
done

#download sample .sra files
for SRA_ID in ${SRA_IDS[@]}; do
    echo "Currently downloading: $SRA_ID"
    prefetch $SRA_ID
done

#convert .sra files into .fastq
#--split-3:
#splits paired-end reads into *_1.fastq and *_2.fastq files, unmated reads are placed into *.fastq
#--skip-technical:
#ignores technical reads (adapters, primers, barcodes, etc)
for SRA_ID in ${SRA_ID[@]}; do
    echo "Generating fastq for: $SRA_ID"
    fasterq-dump $SRA_ID --split-3 --skip-technical --outdir $SRA_ID
    
    #remove the .sra file
    rm -rfv $SRA_ID/*.sra
done

