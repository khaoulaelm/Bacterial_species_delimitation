#!/bin/bash

# Step 1: Extract 16S rRNA genes from Prokka .ffn

mkdir -p 16S_seq
for ffn in prokka2/*/*.ffn; do
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
cat 16S_all/*.fasta > 16S_all_final.fasta

# Step 4: Align with MAFFT
mafft --auto 16S_all_final.fasta > 16S_all_aligned.fasta
