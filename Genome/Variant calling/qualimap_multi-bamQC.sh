#!/bin/bash

#Alignment quality control - QualiMap (multi-BAM QC)
#---
#This script automates the quality control of read alignment for multiple BAM files.


#Dependencies
#---
#The main argument required is the configuration file describing input data (-d).
#This has to be a 2- or 3-column tab-delimted file. The first column should contain
#the sample name and the second column should contain either the absolute or relative
#path to the results of BAM QC analysis or path to the BAM file (if -r mode is activated).
#Additionally the third optional column can provide the condition of the sample. This is an
#optional column. However, if conditions are available they should be provided for each sample.


#define global variables
DIR="/path/to/directory" #working directory

#create required directories
mkdir -p "$DIR/outputs/quality_control/qualimap"
mkdir -p "$DIR/logs/quality_control"


#Multi-sample alignment quality control
#---
#define the config file
config_file="${DIR}/outputs/quality_control/qualimap/multi-bamQC_config.txt"

#define output directory and file
output_dir="${DIR}/outputs/quality_control/qualimap"

#record start time
start_time=$(date +%s)

#run the alignment quality control
qualimap multi-bamqc -d $config_file -outdir $output_dir -r --java-mem-size=8G \
2> "${DIR}/logs/quality_control/multi-bamQC.log"

#record the end time
end_time=$(date +%s)

#store the sample runtime
run_time=$(( end_time - start_time))
minutes=$(( run_time / 60 ))
echo "The multi-bam quality control of read alignments took ${minutes} minutes"
