#!/bin/bash

#Author: Khaoula El Mchachti
# Description: Downloads genome FASTA files from NCBI using assembly accessions.
#Input: accessions.txt (One NCBI assembly accession per line). 
# Output: ncbi_genomes.zip 
# Date: 2026-08-18

echo "===== Genome download started ====="

# Download genomes listed in accessions.txt
datasets download genome accession \
    --inputfile accessions.txt \
    --include genome \
    --filename ncbi_genomes.zip

if [ $? -ne 0 ]; then
    echo "ERROR: NCBI genome download failed."
    exit 1
fi

echo "NCBI genome download completed."
