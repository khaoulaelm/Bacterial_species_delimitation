#!/bin/bash

echo "Starting tree creation from core gene alignment"

# Create Newick tree from core gene alignment
echo "Building Newick tree with FastTree..."
FastTree -nt roary_results/core_gene_alignment.aln > roary_results/tree.newick


