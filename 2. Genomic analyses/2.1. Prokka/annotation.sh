#!/bin/bash

mkdir -p prokka_results

for file in ~/Bacterial_species_delimitation/1. Bacterial strains/genomes/*.fasta; do
    base=$(basename "$file" .fasta)
    prokka --outdir prokka_results/$base --prefix $base --cpus 4 "$file"
done
