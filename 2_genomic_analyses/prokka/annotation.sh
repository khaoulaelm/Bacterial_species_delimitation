#!/bin/bash
# Author: Khaoula El Mchachti
# Project: Bacterial species delimitation
# Description: Annotate genomes using Prokka
#Input: *.fasta genomes
#Output: prokka_results/<strain_name>/ (annotation files)
#Date: 2026-03-02

mkdir -p prokka_results

GENOMES_DIR="$HOME/Bacterial_species_delimitation/1. Bacterial strains/genomes"

for file in "$GENOMES_DIR"/*.fasta; do
    base="$(basename "$file" .fasta)"
    prokka --outdir "prokka_results/$base" --prefix "$base" --cpus 4 "$file"
done

