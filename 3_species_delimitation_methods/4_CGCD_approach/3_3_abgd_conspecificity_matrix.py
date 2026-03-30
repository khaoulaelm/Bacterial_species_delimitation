#!/usr/bin/env python3

"""
Author: Khaoula El Mchachti
Description: Generate the conspecificity matrix by summing ABGD per-gene partition matrices.
Input: ABGD_partition_matrices/ 
Output: ABGD_conspecificity_matrix.csv (pairwise counts of how many genes place two strains in the same group)
Date: 2026-03-20
"""

import os
import pandas as pd
import numpy as np

# Directory where partition matrices are stored
partition_dir = os.path.expanduser("~/Bacterial_species_delimitation/3_species_delimitation_methods/4_CGCD_approach/ABGD_partition_matrices")

# Directory where the file will be saved
output_dir = os.path.expanduser("~/Bacterial_species_delimitation/3_species_delimitation_methods/4_CGCD_approach/ABGD_conspecificity_matrix")

# Create the directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# File path inside the directory
output_path = os.path.join(output_dir, "ABGD_conspecificity_matrix.csv")

print("===== Generating conspecificity matrix =====")


# Get all partition matrix files
files = [f for f in os.listdir(partition_dir) if f.endswith(".csv")]

# Check if partition matrices are available
if not files:
    print(" ERROR: No partition matrices found.")
    exit()

# Use the first file to get all strains
first_matrix = pd.read_csv(os.path.join(partition_dir, files[0]), index_col=0)
strains = list(first_matrix.index)
matrix = pd.DataFrame(0, index=strains, columns=strains, dtype=int)

# Accumulate all matrices
for f in files:
    path = os.path.join(partition_dir, f)
    df = pd.read_csv(path, index_col=0)
    if list(df.index) != strains:
        print(f" WARNING: Strain mismatch in {f}. Skipping.")
        continue
    matrix += df

# Save the final matrix
matrix.to_csv(output_path)
print(f" Conspecificity matrix saved to:\n{output_path}")
