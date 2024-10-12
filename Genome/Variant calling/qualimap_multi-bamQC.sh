#!/bin/bash

#Alignment quality control - QualiMap (multi-BAM QC)
#---
#This script automates the alignment of reads to a reference genome using BWA-MEM.

#define global variables
DIR="/path/to/directory" #working directory

#create required directories
mkdir -p "$DIR/outputs/quality_control/qualimap"
mkdir -p "$DIR/logs/quality_control"


#Multi-BAM alignment quality control
#---
#define the config file
config_file="${DIR}/outputs/quality_control/qualimap/multi-bamQC_config.txt"

#define output directory
output_dir="${DIR}/outputs/quality_control/qualimap"

#record start time
start_time=$(date +%s)

#run the alignment quality control
qualimap multi-bamqc -d $config_file -outdir $output_dir -r \
2> "${DIR}/logs/quality_control/multi-bamQC.log"

#record the end time
end_time=$(date +%s)

#store the sample runtime
run_time=$(( end_time - start_time))
minutes=$(( run_time / 60 ))
echo "The multi-bam quality control of read alignments took ${minutes} minutes"
