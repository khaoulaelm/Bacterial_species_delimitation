#!/bin/bash
# Author: Khaoula El Mchachti
# Description: Annotate genomes using Prokka
#Input: *.fasta genomes
#Output: prokka_results/<strain_name>/ (annotation files)
#Date: 2026-04-03

echo "===== Starting Prokka analysis ====="

mkdir -p prokka_results

GENOMES_DIR="$HOME/Bacterial_species_delimitation/4_large_scale_genome_dataset/CGCD_approach/fasta/renamed"

for file in "$GENOMES_DIR"/*.fna; do
    base=$(basename "$file" .fna)
    prokka --outdir prokka_results/$base --prefix $base --cpus 4 "$file"
done
