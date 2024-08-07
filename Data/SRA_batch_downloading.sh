#!/bin/bash

#Batch downloading SRA files
#---
#This script automates the batch downloading of SRA files using the SRA toolkit.
#SRA-toolkits documentation: https://github.com/ncbi/sra-tools/wiki


#define global variables
DIR="/path/to/directory" #working directory

#define the range of SRA sample id's
SRA_MIN=6298258
SRA_MAX=6298307

#create required directories
mkdir -p "$DIR/samples"


#Download SRA files
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
