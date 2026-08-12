#!/bin/bash
# Author: Khaoula El Mchachti
# Description: Run Roary pan-genome analysis
# Input: Prokka .gff files
# Output: roary_results/
# Date: 2026-03-02
# Last modified: 2026-08-12

echo "===== Starting Roary analysis ====="

# Find the directory
SCRIPT_DIR="$(cd "$(dirname "$BASH_SOURCE[0]}")" && pwd)"

# Project root directory
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Prokka results directory
PROKKA_DIR="$PROJECT_DIR/2_genomic_analyses/1_prokka/prokka_results"

# GFF directory
GFF_DIR="$PROJECT_DIR/2_genomic_analyses/2_roary/gff_files"

# Roary results directory
ROARY_DIR="$PROJECT_DIR/2_genomic_analyses/2_roary/roary_results"

# Prepare the GFF files
echo "Creating GFF directory..."
mkdir -p "$GFF_DIR"
cp "$PROKKA_DIR"/*/*.gff "$GFF_DIR/"

# Run Roary with default parameters
echo "Running Roary (this may take a while)..."
roary -e --mafft -p 8 -f "$ROARY_DIR" "$GFF_DIR"/*.gff
