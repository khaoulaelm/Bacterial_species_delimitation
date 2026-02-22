#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np

# Define the directory where partition matrices are stored
partition_dir = os.path.expanduser("~/partition_matrices")

# Define the output path where the conspecificity matrix will be saved
output_path = os.path.expanduser("~/conspecificity_matrix.csv")

# Get all partition matrix files
files = [f for f in os.listdir(partition_dir) if f.endswith(".csv")]

# Check if partition matrices are available
if not files:
    print(" No partition matrices found.")
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
        print(f" Strain mismatch in {f}. Skipping.")
        continue
    matrix += df

# Save the final matrix
matrix.to_csv(output_path)
print(f" Conspecificity matrix saved to:\n{output_path}")
