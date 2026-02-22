import os
import pandas as pd

# Define directories
base_dir = os.path.expanduser("~/")
partition_dir = os.path.join(base_dir, "asap_partition_matrices")  # ASAP partition matrices
output_path = os.path.join(base_dir, "asap_conspecificity_matrix.csv")

# Get all CSV files in the partition directory
files = [f for f in os.listdir(partition_dir) if f.endswith(".csv")]

if not files:
    print(" No partition matrices found.")
    exit()

# Load the first matrix to initialize structure
first_path = os.path.join(partition_dir, files[0])
first_matrix = pd.read_csv(first_path, index_col=0)
strains = list(first_matrix.index)
matrix = pd.DataFrame(0, index=strains, columns=strains, dtype=int)

# Loop through all partition matrices
for f in files:
    path = os.path.join(partition_dir, f)
    try:
        df = pd.read_csv(path, index_col=0)

        # Check for missing rows/columns
        missing_rows = set(strains) - set(df.index)
        missing_cols = set(strains) - set(df.columns)

        if missing_rows or missing_cols:
            print(f" Skipping {f}:")
            if missing_rows:
                print(f" Missing strains in rows: {sorted(missing_rows)}")
            if missing_cols:
                print(f" Missing strains in columns: {sorted(missing_cols)}")
            continue

        # Reorder rows and columns to match reference
        df = df.loc[strains, strains]

        # Accumulate values
        matrix += df

    except Exception as e:
        print(f" Error processing {f}: {e}")

# Save final conspecificity matrix
matrix.to_csv(output_path)
print(f"\n ASAP conspecificity matrix saved to:\n{output_path}")
