#!/usr/bin/env python3

"""
Author: Khaoula El Mchachti
Description: Align all core gene FASTA files using MAFFT
Input: core_genes_fasta/*.fasta
Output: core_genes_aligned/*_aligned.fasta
Date: 2026-03-02
"""

import os
import subprocess

# Define input/output
input_dir = os.path.expanduser("~/Bacterial_species_delimitation/3_species_delimitation_methods/4_CGCD_approach/core_genes_fasta")
output_dir = os.path.expanduser("~/Bacterial_species_delimitation/3_species_delimitation_methods/4_CGCD_approach/core_genes_aligned")

os.makedirs(output_dir, exist_ok=True)

print("===== Starting MAFFT alignments =====")

for fasta_file in os.listdir(input_dir):
    if fasta_file.endswith(".fasta"):
        input_path = os.path.join(input_dir, fasta_file)
        output_path = os.path.join(output_dir, fasta_file.replace(".fasta", "_aligned.fasta"))

        print(f"Aligning {fasta_file}...")
        subprocess.run(["mafft", "--auto", input_path], stdout=open(output_path, "w"))

print("All alignments saved to:", output_dir)
