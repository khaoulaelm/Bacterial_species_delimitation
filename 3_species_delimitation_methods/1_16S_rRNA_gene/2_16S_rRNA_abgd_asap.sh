#!/bin/bash
# Author: Khaoula El Mchachti
# Description: Run ABGD and ASAP species delimitation analyses from aligned 16S sequences
# Input: 16S_aligned.fasta
# Output: abgd_output/, asap_output/
# Date: 2026-03-02

#ABGD
echo "Running ABGD species delimitation..."

mkdir -p abgd_output

~/Bacterial_species_delimitation/3_species_delimitation_methods/ABGD/abgd \
  -a \
  -o abgd_output/ \
  -d JC69 \
  ~/Bacterial_species_delimitation/3_species_delimitation_methods/1_16S_rRNA_gene/16S_aligned.fasta

echo "ABGD finished"


# ASAP

echo "Running ASAP species delimitation..."

mkdir -p asap_output

~/Bacterial_species_delimitation/3_species_delimitation_methods/ASAP/asap \
  -u\
  -o asap_output/ \
  ~/Bacterial_species_delimitation/3_species_delimitation_methods/1_16S_rRNA_gene/16S_aligned.fasta
  
echo "ASAP finished"
