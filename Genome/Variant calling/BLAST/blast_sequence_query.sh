#!/bin/bash

#BLAST search
#---
#Script to automate BLAST searches to compare query nucleotide sequences against a reference genome.


#Table of contents:
#---
#1. Create organisms reference sequence database
#2. BLAST transcript sequence against reference organism


#define global variables
DIR="/path/to/directory" #working directory
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.4_UCB_Xtro_10.0/GCF_000004195.4_UCB_Xtro_10.0_genomic.fna.gz" #reference sequence FTP link
REF_ORG="UCB_Xtro_10.0_db" #reference organism name
TRANSCRIPT_IDS=("NM_204089.1" "NM_204090.1" "NM_204091.1") #transcript ids (add more if needed)

#create required directories
mkdir -p "${DIR}/transcripts"
mkdir -p "${DIR}/outputs"


#1. Create organisms reference sequence database (if not already created)
#---
if false; then
	#download the reference sequence
	echo "Downloading reference sequence..."
	wget -P $DIR $REF_URL

	#define the compressed reference sequence file name
	ref_gz=$(find "${DIR}" -type f -name "*.fna.gz" | head -n 1)

	#uncompress the file
	echo "Unzipping reference sequence..."
	gunzip -d $ref_gz

	#define the reference sequence file name
	ref=$(find "$DIR" -type f -name "*.fna" | head -n 1)

	#define the BLAST database name
	blast_db="${DIR}/${REF_ORG}"

	#create the BLAST database	
	echo "Creating BLAST database for reference sequence..."
	makeblastdb -in $ref -dbtype nucl -out $blast_db
	echo "BLAST reference sequence successfully created!"
fi


#2. BLAST transcript sequence against reference organism
#---
#define function for parrallelization
run_blast() {
	#define local transcript variables
	local transcript_id=$1
	local transcript_fasta="${DIR}/transcripts/${transcript_id}.fa"
	local transcript_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${transcript_id}&rettype=fasta&retmode=text"

	#define the local BLAST variables
	local blast_db="${DIR}/${REF_ORG}"
	local blast_result="${DIR}/outputs/blastn_${transcript_id}.txt"

	#download the transcript sequence
	wget -O $transcript_fasta $transcript_url

	#run BLAST if the transcript download was successful
    if [ -f $transcript_fasta ]; then
        blastn -query $transcript_fasta -db $blast_db -out $blast_result
    else
        echo "Error: Failed to download ${transcript_id} sequence"
    fi
}

#export function and variables for GNU parallel
export -f run_blast
export DIR REF_ORG

#run GNU parallel
parallel -j 8 run_blast ::: "${TRANSCRIPT_IDS[@]}"
