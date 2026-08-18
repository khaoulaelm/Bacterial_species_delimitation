#!/usr/bin/env python3

"""
Author: Khaoula El Mchachti
Description: Extract all core genes (present in all strains)
Input: roary_results/gene_presence_absence.csv, prokka_results/
Output: /core_genes_fasta/<gene_name>.fasta, missing_genes_log.csv
Date: 2026-03-20
Last modified: 2026-08-18
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

# Missing-data log
missing_log_path = os.path.join(
    output_dir,
    "missing_genes_log.csv"
)    

# Make sure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load roary matrix
df = pd.read_csv(csv_path, low_memory=False)

# Get strain names from column headers (strain columns start from column 15)
strain_columns = df.columns[14:]  
number_of_strains = len(strain_columns)
print(f"Number of isolates detected: {number_of_strains}")

# Filter core genes (present in all strains)
core_genes = df[df["No. isolates"].astype(float) == number_of_strains]

print(f" Found {len(core_genes)} strict core genes (present in all {len(strain_columns)} strains).")

# Initialize missing data logger 
missing_data = []

# Extract core genes sequences
for idx, row in core_genes.iterrows():
    gene_name = row["Gene"] if pd.notna(row["Gene"]) else f"group_{idx+1}"
    output_fasta = os.path.join(output_dir, f"{gene_name}.fasta")
    sequences_written = 0
    
    # Create FASTA file for this core gene
    with open(output_fasta, "w") as out_fasta:
        for strain in strain_columns:
            gene_ids_str = row[strain]
            
            # Check whether Roary contains a gene ID
            if pd.isna(gene_ids_str):
                missing_data.append((gene_name, strain, "Missing gene ID in Roary"))
                continue

            # Split IDs by semicolon or spaces
            gene_ids = [g.strip() for g in str(gene_ids_str).replace(";", " ").split() if g.strip()]
            
            # Path to Prokka .ffn file
            strain_ffn = os.path.join(prokka_base_dir, strain, f"{strain}.ffn")
            
            # Check .ffn file
            if not os.path.isfile(strain_ffn):
                missing_data.append((gene_name, strain, "Missing .ffn file"))
                continue

            # Read Prokka nucleotide sequences
            records = list(SeqIO.parse(strain_ffn, "fasta"))
            matched = False
            
            # Find the corresponding gene sequence
            for gene_id in gene_ids:
                for record in records:
                    if gene_id in record.id:  # Match ONLY by ID
                        record.id = strain
                        record.description = ""
                        SeqIO.write(record, out_fasta, "fasta")
                        print(f"{gene_name}  {strain} (selected {record.id})")
                        sequences_written += 1
                        matched = True
                        break
                if matched:
                    break
            
            # Gene ID was not found in Prokka .ffn
            if not matched:
                missing_data.append((gene_name, strain, "Gene ID not found in .ffn"))

    # Remove empty FASTA files
    if sequences_written == 0:
        os.remove(output_fasta)

# Save missing gene log 
missing_df = pd.DataFrame(missing_data, columns=["Gene", "Strain", "Reason"])
missing_df.to_csv(missing_log_path, index=False)

print(f"\n Extraction complete. FASTA files saved in: {output_dir}")
print(f"Missing gene log saved to: {missing_log_path}")
