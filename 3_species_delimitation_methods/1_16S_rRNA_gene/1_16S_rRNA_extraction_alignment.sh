#!/bin/bash
# Author: Khaoula El Mchachti
# Description: Extract 16S rRNA sequences from Prokka annotations and perform MAFFT alignment
# Input: Prokka .ffn files
# Output: 16S_sequences.fasta, 16S_aligned.fasta
# Date: 2026-03-02
#Last modified: 2026-08-12

echo "===== Extracting and aligning 16S sequences ====="

# Find the directory
SCRIPT_DIR="$(cd "$(dirname "$BASH_SOURCE[0]}")" && pwd)"

#Project root directory
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Prokka results directory
PROKKA_DIR="$PROJECT_DIR/2_genomic_analyses/1_prokka/prokka_results"

# 16S output directory
OUTPUT_DIR="$PROJECT_DIR/3_species_delimitation_methods/1_16S_rRNA_gene"

# Step 1: Extract 16S rRNA genes from Prokka .ffn
echo "Extracting 16S sequences from .ffn files..."
mkdir -p "$OUTPUT_DIR/16S_seq"

for ffn in "$PROKKA_DIR"/*/*.ffn; do
    base=$(basename "$ffn" .ffn)
    seqkit grep -nrp "16S.*ribosomal.*RNA" "$ffn" > "16S_seq/${base}.fasta"
done

# Step 2: Rename headers to keep all copies traceable: >strain_copy1, >strain_copy2, ...
mkdir -p "$OUTPUT_DIR/16S_all"
for file in "$OUTPUT_DIR/16S_seq"/*.fasta; do
    base=$(basename "$file" .fasta)  # strain name
    
    awk -v strain="$base" '
        /^>/ {++i; print ">"strain"_copy"i}
        !/^>/ {print}
    ' "$file" > "16S_all/${base}.fasta"
done

# Step 3: Concatenate all sequences
cat "$OUTPUT_DIR/16S_all"/*.fasta > "$OUTPUT_DIR/16S_sequences.fasta"

# Step 4: Align with MAFFT
echo "Running MAFFT alignment..."
mafft --auto "$OUTPUT_DIR/16S_sequences.fasta" > "$OUTPUT_DIR/16S_aligned.fasta"

# Step 5: Clean intermediate files (keep only final outputs)
echo "Cleaning intermediate files..."
rm -rf "$OUTPUT_DIR/16S_seq" "$OUTPUT_DIR/16S_all"

echo "Done"
echo "Final files:"
echo "16S_sequences.fasta"
echo "16S_aligned.fasta"
