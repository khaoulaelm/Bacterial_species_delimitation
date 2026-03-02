#!/bin/bash

echo "Starting Roary plotting"

# Plot Roary results
echo "Plotting Roary results..."
python3 roary_plots.py \
  --labels \
  roary_results/tree.newick \
  roary_results/gene_presence_absence.csv
