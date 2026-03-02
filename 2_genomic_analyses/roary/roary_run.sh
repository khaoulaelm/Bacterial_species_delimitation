#!/bin/bash

# Prepare the GFF files
echo "Creating GFF directory..."
mkdir -p gff_files
cp ~/Bacterial_species_delimitation/2_genomic_analyses/prokka/prokka_results/*/*.gff gff_files/

# Run roary 
# Create output directory

echo "Running Roary (this may take a while)..."


# Run Roary with default parameters
roary -e --mafft -p 8 -f roary_results gff_files/*.gff

