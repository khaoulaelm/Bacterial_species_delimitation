#!/bin/bash
# Author: Khaoula El Mchachti
# Project: Bacterial species delimitation
# Description: Annotate genomes using Prokka
#Input: *.fasta genomes
#Output: prokka_results/<strain_name>/ (annotation files)
#Date: 2026-03-02
#Last modified: 2026-08-12

echo "===== Starting Prokka analysis ====="

# Find the directory
SCRIPT_DIR="$(cd "$(dirname "$BASH_SOURCE[0]}")" && pwd)"

#Project root directory
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"


GENOMES_DIR="$PROJECT_DIR/1_bacterial_strains/genomes"

mkdir -p prokka_results

for file in "$GENOMES_DIR"/*.fasta; do
    base="$(basename "$file" .fasta)"
    prokka --outdir "prokka_results/$base" --prefix "$base" --cpus 4 "$file"
done
