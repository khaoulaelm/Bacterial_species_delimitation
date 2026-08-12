#!/bin/bash
# Author: Khaoula El Mchachti
# Description: Extract 16S rRNA sequences from Prokka annotations and perform MAFFT alignment
# Input: Prokka .ffn files
# Output: 16S_sequences.fasta, 16S_aligned.fasta
# Date: 2026-03-02

echo "===== Extracting and aligning 16S sequences ====="

# Step 1: Extract 16S rRNA genes from Prokka .ffn
echo "Extracting 16S sequences from .ffn files..."
mkdir -p 16S_seq

for ffn in ~/Bacterial_species_delimitation/2_genomic_analyses/prokka/prokka_results/*/*.ffn; do
    base=$(basename "$ffn" .ffn)
    seqkit grep -nrp "16S.*ribosomal.*RNA" "$ffn" > "16S_seq/${base}.fasta"
done

# Step 2: Rename headers to keep all copies traceable: >strain_copy1, >strain_copy2, ...
mkdir -p 16S_all
for file in 16S_seq/*.fasta; do
    base=$(basename "$file" .fasta)  # strain name
    awk -v strain="$base" '
        /^>/ {++i; print ">"strain"_copy"i}
        !/^>/ {print}
    ' "$file" > "16S_all/${base}.fasta"
done

# Step 3: Concatenate all sequences
cat 16S_all/*.fasta > 16S_sequences.fasta

# Step 4: Align with MAFFT
echo "Running MAFFT alignment..."
mafft --auto 16S_sequences.fasta > 16S_aligned.fasta

# Step 5: Clean intermediate files (keep only final outputs)
echo "Cleaning intermediate files..."
rm -rf 16S_seq 16S_all

echo "Done"
echo "Final files:"
echo "  16S_sequences.fasta"
echo "  16S_aligned.fasta"
