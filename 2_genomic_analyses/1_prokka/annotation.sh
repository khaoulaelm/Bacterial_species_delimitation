#!/bin/bash
# Author: Khaoula El Mchachti
# Project: Bacterial species delimitation
# Description: Annotate genomes using Prokka
#Input: *.fna genomes
#Output: prokka_results/<strain_name>/ (annotation files)
#Date: 2026-03-02
#Last modified: 2026-08-12

echo "===== Starting Prokka analysis ====="

# Find the directory
SCRIPT_DIR="$(cd "$(dirname "$BASH_SOURCE[0]}")" && pwd)"

# Project root directory
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Genome directory
GENOMES_DIR="$PROJECT_DIR/1_bacterial_strains/genomes"

# Prokka output directory
PROKKA_DIR="$SCRIPT_DIR/prokka_results"

mkdir -p "$PROKKA_DIR"

for file in "$GENOMES_DIR"/*.fna; do
    base="$(basename "$file" .fna)"
    prokka --outdir "$PROKKA_DIR/$base" --prefix "$base" --cpus 4 "$file"
done

echo "Prokka analysis completed"
echo "Output drectory: $PROKKA_DIR"
