#!/bin/bash
# Author: Khaoula El Mchachti
# Description: Run FastANI all-vs-all genome comparison
# Input: genomes_paths_list.txt 
# Output: fastani_results.csv
# Date: 2026-03-02

echo "===== Running FastANI ====="

mkdir -p ani_results

GENOMES_DIR="${1:-$HOME/Bacterial_species_delimitation/1_bacterial_strains/genomes}"
GENOME_LIST="ani_results/genomes_paths_list.txt"

echo "Generating genome list automatically..."
echo "NOTE: genomes_paths_list.txt is created dynamically."
echo "Do NOT edit it manually. Absolute paths are required for FastANI."

# Create genome list with absolute paths
find "$GENOMES_DIR" -type f -name "*.fasta" > "$GENOME_LIST"

echo "Genome list created at: $GENOME_LIST"

# Run FastANI
echo "Running FastANI..."
fastANI \
  --ql "$GENOME_LIST" \
  --rl "$GENOME_LIST" \
  -o ani_results/fastani_results.csv

echo "Output: ani_results/fastani_results.csv"
