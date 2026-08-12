#!/bin/bash
# Author: Khaoula El Mchachti
# Description: Plotting Roary results using roary_plots.py
# Input: tree.newick, gene_presence_absence.csv
# Output: Roary plot files (plots)
# Date: 2026-03-02

echo "===== Starting Roary plotting ======"

# Plot Roary results
echo "Plotting Roary results..."
python3 roary_plots.py \
  --labels \
  roary_results/tree.newick \
  roary_results/gene_presence_absence.csv
