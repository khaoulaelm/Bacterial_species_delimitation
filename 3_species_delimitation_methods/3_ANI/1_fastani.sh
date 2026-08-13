#!/bin/bash
# Author: Khaoula El Mchachti
# Description: Run FastANI all-vs-all genome comparison
# Input: genomes_paths_list.txt 
# Output: fastani_results.csv
# Date: 2026-03-02
# Last modified: 2026-08-13

# Find the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "$BASH_SOURCE[0]}")" && pwd)"

# Project root directory
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Genome directory
GENOMES_DIR="$PROJECT_DIR/1_bacterial_strains/genomes"

# ANI results directory
ANI_DIR="$SCRIPT_DIR/ani_results"

# Genome list
GENOME_LIST="$ANI_DIR/genomes_paths_list.txt"


echo "===== Running FastANI ====="

echo "Please activate the conda environment containing FastANI before runing this script."
mkdir -p "$ANI_DIR"


echo "Generating genome list automatically..."

# Create genome list with absolute paths
find "$GENOMES_DIR" -type f -name "*.fasta" > "$GENOME_LIST"

echo "Genome list created at: $GENOME_LIST"

# Run FastANI
echo "Running FastANI..."
fastANI \
  --ql "$GENOME_LIST" \
  --rl "$GENOME_LIST" \
  -o ani_results/fastani_results.csv
  
echo "FastANI finished"
echo "Output: $ANI_DIR/fastani_results.csv"
