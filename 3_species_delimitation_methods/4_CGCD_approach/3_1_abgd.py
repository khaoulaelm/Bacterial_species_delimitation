#!/usr/bin/env python3

"""
Author: Khaoula El Mchachti
Description: Run ABGD on each aligned core-gene (one gene at a time)
Input: core_genes_aligned/
Output: ABGD_results/<gene_name>/ (ABGD results per gene), failed_genes.txt (log file listing failures)
Date: 2026-03-20
Last modified: 2026-08-16
"""

import os
import subprocess

# Find the directory containing this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Project root directory
PROJECT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "../.."))

# Path to the ABGD executable
abgd_exec = os.path.join(
    PROJECT_DIR,
    "3_species_delimitation_methods",
    "ABGD",
    "abgd"
)  

# Directory containing the aligned core genes (in FASTA format)
input_dir = os.path.join(
    PROJECT_DIR,
    "3_species_delimitation_methods",
    "4_CGCD_approach",
    "core_genes_aligned"
)

# Directory where the ABGD results will be stored
output_dir = os.path.join(
    PROJECT_DIR,
    "3_species_delimitation_methods", 
    "4_CGCD_approach",   
    "ABGD_results"
)

# ABGD is an older executable and may have problems withh long absolute paths.
# Therefore, I convert the automatically detected paths to relative paths before passing them to ABGD.

abgd_exec_relative = os.path.relpath(
    abgd_exec,
    SCRIPT_DIR
) 

# Ensure output root exists
os.makedirs(output_dir, exist_ok=True)

failed_genes_file = os.path.join(output_dir, "failed_genes.txt")

print("===== Running ABGD on core genes =====")
print("Please make sure that the appropriate environment is activated.")

# Make ABGD executable 
subprocess.run(["chmod", "+x", abgd_exec])

# Open the failed genes log file in write mode
with open(failed_genes_file, 'w') as failed_genes_log:
    failed_genes_log.write("Failed genes:\n")  # Write a header to the log file

    # Loop over each file in the input directory
    for fname in os.listdir(input_dir):
        # Skip files that don't end with .fasta (we only want the FASTA files)
        if not fname.endswith(".fasta"):
            continue

        # Get the full path of the current gene alignment file
        gene_path = os.path.join(input_dir, fname)
        
        # Extract the gene name     
        gene_name = os.path.splitext(fname)[0]

        # Create a unique folder for each gene output
        gene_output_dir = os.path.join(output_dir, gene_name)
        os.makedirs(gene_output_dir, exist_ok=True)
        
        # ABGD is an older executable and may have problems withh long absolute paths.
        # Therefore, I convert the automatically detected paths to relative paths before passing them to ABGD.
        gene_path_relative = os.path.relpath(
            gene_path,
            SCRIPT_DIR
        ) 
                
        gene_output_relative = os.path.relpath(
            gene_output_dir,
            SCRIPT_DIR
        ) 
        # Build the ABGD command
        cmd = [
            abgd_exec_relative,
            "-a",                           # Output all files
            "-o", gene_output_relative,     # Output to this gene-specific folder
            gene_path_relative              # Input gene alignment
        ]
        
        # Print a message indicating which gene ABGD is currently processing
        print(f"Running ABGD on: {gene_name}")
        
        # Run the ABGD command and wait for it to complete   
        print("Command being executed:")
        print(" ".join(cmd)) 
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=SCRIPT_DIR)
           
        # Check the ABGD output to see if it completed successfully
        if result.returncode != 0 or "ERROR" in result.stderr:
            # If ABGD fails (non-zero exit code or ERROR in stderr), log the failure
            print(f"ABGD failed for {gene_name}. Adding to failed list.")
            print(result.stderr)

            failed_genes_log.write(f"{gene_name}\n")  # Log the failed gene name
        else:
            print(f"ABGD completed successfully for {gene_name}.")

# After all genes are processed, print a final message
print("All genes processed with ABGD. Failed genes are logged in:", failed_genes_file)
