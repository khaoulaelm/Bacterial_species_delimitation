#!/bin/bash
# Author: Khaoula El Mchachti
# Description: Generate phylogenetic tree from core gene alignment using FastTree 
# Input: core_gene_alignment.aln
# Output: tree.newick
# Date: 2026-03-02

echo "===== Starting tree creation from core gene alignment ====="

# Create Newick tree from core gene alignment
echo "Building Newick tree with FastTree..."
FastTree -nt roary_results/core_gene_alignment.aln > roary_results/tree.newick


