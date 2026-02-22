#!/usr/bin/env python3
import os
import pandas as pd
from collections import defaultdict

# Define directories for ABGD results and where the partition matrices will be saved
abgd_dir = os.path.expanduser("~/ABGD_results")
output_dir = os.path.expanduser("~/partition_matrices")

# Make sure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Loop through each gene folder in the ABGD results directory
for gene_folder in os.listdir(abgd_dir):
    gene_path = os.path.join(abgd_dir, gene_folder)
    if not os.path.isdir(gene_path):
        continue

    # Look for .res.cvs file
    res_file = [f for f in os.listdir(gene_path) if f.endswith(".res.cvs")]
    
    # If no .res.cvs file is found, skip the current folder and print a message
    if not res_file:
        print(f" No .res.cvs file for {gene_folder}")
        continue

    # Read the .res.cvs file 
    res_df = pd.read_csv(os.path.join(gene_path, res_file[0]), sep="\t")
    if "nbSubsetRecursive" not in res_df.columns:
        print(f" No 'nbSubsetRecursive' column in {res_file[0]}")
        continue

    # Most frequent recursive value
    freq_counts = res_df["nbSubsetRecursive"].value_counts()
    top_recursive = freq_counts.index[0]

    # Get last partition index with that value
    matches = res_df[res_df["nbSubsetRecursive"] == top_recursive]
    partition_index = matches.index[-1] + 1  # part.# starts from 1
    
    # Construct the partition file path based on the partition index
    part_file = os.path.join(gene_path, f"{gene_folder}.part.{partition_index}.txt")
    
    # If the partition file is not found, skip to the next gene
    if not os.path.isfile(part_file):
        print(f" Partition file not found: {part_file}")
        continue

    # Parse group file to get groupings
    groups = defaultdict(list)
    with open(part_file) as f:
        for line in f:
            if "Group[" in line:
                parts = line.strip().split("id: ")
                if len(parts) > 1:
                    strains = parts[1].split()
                    group_id = line.split("Group[")[1].split("]")[0].strip()
                    groups[group_id].extend(strains)

    # Get all strains
    all_strains = sorted({strain for group in groups.values() for strain in group})

    # Build partition matrix
    matrix = pd.DataFrame(0, index=all_strains, columns=all_strains, dtype=int)
    for group in groups.values():
        for i in group:
            for j in group:
                matrix.loc[i, j] = 1

    # Save the matrix to a CSV file in the output directory
    matrix.to_csv(os.path.join(output_dir, f"{gene_folder}.csv"))
    
    # Print a message indicating the matrix has been saved for the current gene
    print(f" Saved matrix for {gene_folder}")

# Print a message indicating the matrix has been saved for the current gene
print(" All ABGD partition matrices generated.")
