#!/bin/bash
# Author: Khaoula El Mchachti
# Description: Run ABGD and ASAP species delimitation analyses from aligned core-genome
# Input: core_gene_alignment.aln
# Output: abgd_output/, asap_output/
# Date: 2026-03-02
# Last modified: 2026-08-12

# Find the directory
SCRIPT_DIR="$(cd "$(dirname "$BASH_SOURCE[0]}")" && pwd)"

# Project root directory
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# ABGD
ABGD_DIR="$PROJECT_DIR/3_species_delimitation_methods/ABGD/abgd"

# Core genome alignment
ALIGNMENT="$PROJECT_DIR/2_genomic_analyses/2_roary/roary_results/core_gene_alignment.aln"

# ASAP
ASAP_DIR="$PROJECT_DIR/3_species_delimitation_methods/ASAP/asap"

#ABGD
echo "Running ABGD species delimitation..."

# Make ABGD executable 
chmod +x "$ABGD_DIR"

mkdir -p abgd_output

"$ABGD_DIR" \
  -a \
  -o abgd_output \
  -d JC69 \
  "$ALIGNMENT"

echo "ABGD finished"

echo "Running ASAP species delimitation..."

# Make ASAP executable 
chmod +x "$ASAP_DIR"

mkdir -p asap_output

"$ASAP_DIR" \
  -u \
  -o asap_output \
  "$ALIGNMENT"
  
echo "ASAP finished"
