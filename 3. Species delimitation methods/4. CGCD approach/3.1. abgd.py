#!/usr/bin/env python3
import os
import subprocess

# Path to the ABGD executable
abgd_exec = os.path.expanduser("~/Abgd/abgd")  

# Directory containing the aligned core genes (in FASTA format)
input_dir = os.path.expanduser("~/core_genes_aligned")

# Directory where the ABGD results will be stored
main_output_dir = os.path.expanduser("~/ABGD_results")

# Path to the file that will store the list of failed genes
failed_genes_file = os.path.expanduser("~/failed_genes.txt")

# Ensure output root exists
os.makedirs(main_output_dir, exist_ok=True)

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
        gene_output_dir = os.path.join(main_output_dir, gene_name)
        os.makedirs(gene_output_dir, exist_ok=True)

        # Build the ABGD command
        cmd = [
            abgd_exec,
            "-a",                      # Output all files
            "-o", gene_output_dir,     # Output to this gene-specific folder
            gene_path                  # Input gene alignment
        ]
        
        # Print a message indicating which gene ABGD is currently processing
        print(f"Running ABGD on: {gene_name}")
        
        # Run the ABGD command and wait for it to complete    
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check the ABGD output to see if it completed successfully
        if result.returncode != 0 or "ERROR" in result.stderr:
            # If ABGD fails (non-zero exit code or ERROR in stderr), log the failure
            print(f"ABGD failed for {gene_name}. Adding to failed list.")
            failed_genes_log.write(f"{gene_name}\n")  # Log the failed gene name
        else:
            print(f"ABGD completed successfully for {gene_name}.")

# After all genes are processed, print a final message
print("All genes processed with ABGD. Failed genes are logged in:", failed_genes_file)
