import os
import csv
import numpy as np

# Define paths
base_dir = os.path.expanduser("~/ASAP_analysis")
partition_dir = os.path.join(base_dir, "partition_matrices")  # Directory with partition matrices (successful genes)
strains_file = os.path.expanduser("~/strains.txt")  # File with strains
output_dir = os.path.join(base_dir, "conspecificity_matrix")
os.makedirs(output_dir, exist_ok=True)

# Load 47 strain names
with open(strains_file, "r") as f:
    strains = [s.strip() for s in f.read().strip().split(",") if s.strip()]

# Initialize an empty matrix to accumulate conspecificity values
conspecificity_matrix = np.zeros((len(strains), len(strains)))

# Step 1: Process each partition matrix from successful genes
partition_files = [f for f in os.listdir(partition_dir) if f.endswith(".csv")]

for partition_file in partition_files:
    matrix_path = os.path.join(partition_dir, partition_file)
    
    # Read the partition matrix
    with open(matrix_path, "r") as f:
        reader = csv.reader(f)
        matrix = [row[1:] for row in reader][1:]  # Exclude first row and column (strain names)
    
    # Convert to numpy array of integers
    matrix = np.array(matrix, dtype=int)
    
    # Add the matrix to the conspecificity matrix
    conspecificity_matrix += matrix

# Step 2: Save the combined conspecificity matrix to a CSV file
output_path = os.path.join(output_dir, "conspecificity_matrix.csv")
with open(output_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([""] + strains)  # First row with strain names
    for strain, row in zip(strains, conspecificity_matrix):
        writer.writerow([strain] + row.tolist())  # First column with strain names

print(f"Conspecificity matrix saved to: {output_path}")
