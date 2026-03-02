#!/bin/bash
# Author: Khaoula El Mchachti
# Description: Run Roary pan-genome analysis
# Input: Prokka .gff files
# Output: roary_results/
# Date: 2026-03-02

echo "===== Starting Roary analysis ====="

# Prepare the GFF files
echo "Creating GFF directory..."
mkdir -p gff_files
cp ~/Bacterial_species_delimitation/2_genomic_analyses/prokka/prokka_results/*/*.gff gff_files/

# Run Roary with default parameters
echo "Running Roary (this may take a while)..."
roary -e --mafft -p 8 -f roary_results gff_files/*.gff

