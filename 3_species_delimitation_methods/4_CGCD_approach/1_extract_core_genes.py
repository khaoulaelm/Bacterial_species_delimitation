#!/usr/bin/env python3

"""
Author: Khaoula El Mchachti
Description: Extract all core genes (present in all strains)
Input: gene_presence_absence.csv, prokka_results/
Output: /core_genes_fasta/<gene_name>.fasta (one file per core gene; headers renamed to strain)
Date: 2026-03-20
Last modified: 2026-08-16
"""

import os
import pandas as pd
from Bio import SeqIO

print("===== Extracting core gene sequences =====")

# Find the directory containing this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Project root directory
PROJECT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "../.."))

# Paths
csv_path = os.path.join(
    PROJECT_DIR,
    "2_genomic_analyses",
    "2_roary",
    "roary_results",
    "gene_presence_absence.csv"
    )
prokka_base_dir = os.path.join(
    PROJECT_DIR,
    "2_genomic_analyses",
    "1_prokka",
    "prokka_results"
    )
output_dir = os.path.join(
    PROJECT_DIR,
    "3_species_delimitation_methods",
    "4_CGCD_approach",
    "core_genes_fasta"
    )

# Make sure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load CSV
df = pd.read_csv(csv_path)

# Get strain names from column headers
strain_columns = df.columns[14:]  

# Filter core genes (present in all 47 isolates)
core_genes = df[df["No. isolates"] == 47]

print(f"Found {len(core_genes)} core genes")

# Iterate through each core gene
for idx, row in core_genes.iterrows():
    gene_name = row["Gene"] if pd.notna(row["Gene"]) else f"group_{idx+1}"
    output_fasta = os.path.join(output_dir, f"{gene_name}.fasta")

    with open(output_fasta, "w") as out_fasta:
        for strain in strain_columns:
            gene_ids = row[strain]
            if pd.isna(gene_ids):
                continue

            for gene_id in str(gene_ids).split(';'):
                gene_id = gene_id.strip()
                if not gene_id:
                    continue

                strain_ffn_path = os.path.join(prokka_base_dir, strain, f"{strain}.ffn")
                if not os.path.isfile(strain_ffn_path):
                    print(f" Missing file for {strain}: {strain_ffn_path}")
                    continue

                for record in SeqIO.parse(strain_ffn_path, "fasta"):
                    if record.id == gene_id:
                        record.id = strain  # Rename to strain name
                        record.description = ""
                        SeqIO.write(record, out_fasta, "fasta")
                        break

print("Extraction complete. FASTA files saved in:", output_dir)
