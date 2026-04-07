#!/usr/bin/env python3

"""
Author: Khaoula El Mchachti
Description: Run ASAP on each aligned core-gene (one gene at a time)
Input: core_genes_aligned/
Output: ASAP_results/<gene_name>/ (ASAP results per gene), failed_genes.txt (log file listing failures)
Date: 2026-03-20
"""

import os
import subprocess

# Paths
asap_exec = os.path.expanduser("~/Bacterial_species_delimitation/3_species_delimitation_methods/ASAP/asap")
input_dir = os.path.expanduser("~/Bacterial_species_delimitation/3_species_delimitation_methods/4_CGCD_approach/core_genes_aligned")  
output_dir = os.path.expanduser("~/Bacterial_species_delimitation/3_species_delimitation_methods/4_CGCD_approach/ASAP_results")
failed_log = os.path.join(output_dir, "failed_genes.txt")

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Initialize failed log
with open(failed_log, "w") as flog:
    flog.write("Failed Genes:\n")

# Loop through aligned files
for fname in os.listdir(input_dir):
    if not fname.endswith(".fasta"):
        continue

    gene_name = fname.replace(".fasta", "")
    in_path = os.path.join(input_dir, fname)
    out_dir = os.path.join(output_dir, gene_name)
    spart_file = os.path.join(out_dir, f"{fname}.spart")  # Matches ASAP naming

    # Skip already completed genes
    if os.path.exists(spart_file):
        print(f"Skipping {gene_name} (already processed)")
        continue

    os.makedirs(out_dir, exist_ok=True)
    print(f"\nRunning ASAP for {gene_name}...")

    try:
        result = subprocess.run(
            [asap_exec, "-a", "-o", out_dir, in_path],
            capture_output=True,
            text=True,
            timeout=90,  #  Optional timeout safeguard
            cwd=out_dir 
        )

        if result.returncode != 0 or "ASAP     failed" in result.stderr or "invalid pointer" in result.stderr:
            print(f" ASAP failed for {gene_name}")
            with open(failed_log, "a") as flog:
                flog.write(f"{gene_name}\n")
        else:
            print(f" ASAP completed for {gene_name}")

    except subprocess.TimeoutExpired:
        print(f" Timeout: ASAP took too long on {gene_name}. Skipping.")
        with open(failed_log, "a") as flog:
            flog.write(f"{gene_name} (timeout)\n")

    except Exception as e:
        print(f" Exception for {gene_name}: {e}")
        with open(failed_log, "a") as flog:
            flog.write(f"{gene_name} (exception)\n")

print("All genes processed with ABGD. Failed genes are logged in:", failed_genes_file)
