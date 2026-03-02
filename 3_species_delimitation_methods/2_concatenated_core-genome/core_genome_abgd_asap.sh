#!/bin/bash
# Author: Khaoula El Mchachti
# Description: Run ABGD and ASAP species delimitation analyses from aligned core-genome
# Input: core_gene_alignment.aln
# Output: abgd_output/, asap_output/
# Date: 2026-03-02

#ABGD

echo "Running ABGD species delimitation..."

mkdir -p abgd_output

~/Bacterial_species_delimitation/3_species_delimitation_methods/ABGD/abgd \
  -a \
  -o abgd_output/ \
  -d JC69 \
  ~/Bacterial_species_delimitation/2_genomic_analyses/roary/roary_results/core_gene_alignment.aln

echo "ABGD finished"


#ASAP
echo "Running ASAP species delimitation..."

mkdir -p asap_output

~/Bacterial_species_delimitation/3_species_delimitation_methods/ASAP/asap \
  -u\
  -o asap_output/ \
  ~/Bacterial_species_delimitation/2_genomic_analyses/roary/roary_results/core_gene_alignment.aln
  
echo "ASAP finished"
