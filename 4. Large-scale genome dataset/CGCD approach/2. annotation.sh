#!/bin/bash

# Annotation with Prokka

mkdir -p prokka_results

for file in ~renamed/*.fna; do
    base=$(basename "$file" .fna)
    prokka --outdir prokka_results/$base --prefix $base --cpus 4 "$file"
done
