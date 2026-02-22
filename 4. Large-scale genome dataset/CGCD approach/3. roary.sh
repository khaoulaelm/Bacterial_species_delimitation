#!/bin/bash

# Prepare the GFF files
mkdir -p gff_files
cp prokka_results/*/*.gff gff_files/

# Run roary 
# Create output directory
mkdir -p roary_results

# Run Roary with default parameters
roary -e --mafft -p 8 -f roary_results gff_files/*.gff

# Create Newick tree
FastTree -nt roary_results/core_gene_alignment.aln > tree.newick

# Plot results
python3 roary_plots.py --labels roary_results/tree.newick roary_results/gene_presence_absence.csv
