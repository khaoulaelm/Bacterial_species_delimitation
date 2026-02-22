import os
import csv

# Base directory for ASAP results
asap_dir = os.path.expanduser("~/ASAP_results")
output_dir = os.path.join(asap_dir, "partition_matrices")
os.makedirs(output_dir, exist_ok=True)

# Load strain names from strains.txt
with open(os.path.expanduser("~/strains.txt")) as f:
    strains = [s.strip() for s in f.read().split(",") if s.strip()]

# Create a matrix where 1 = same group, 0 = different
def build_partition_matrix(strains, group_dict):
    matrix = []
    for strain1 in strains:
        row = []
        for strain2 in strains:
            same_group = int(group_dict.get(strain1) == group_dict.get(strain2))
            row.append(same_group)
        matrix.append(row)
    return matrix

# Save the matrix to CSV with strain labels
def save_partition_matrix(strains, matrix, output_path):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([""] + strains)
        for strain, row in zip(strains, matrix):
            writer.writerow([strain] + row)

# Get full list of gene directories
gene_dirs = [d for d in os.listdir(asap_dir) if os.path.isdir(os.path.join(asap_dir, d))]

# Process each gene
total = len(gene_dirs)
for gene in gene_dirs:
    res_file = os.path.join(asap_dir, gene, f"{gene}.fasta.res.cvs")
    if not os.path.exists(res_file):
        continue

    try:
        with open(res_file) as f:
            lines = f.readlines()
        if len(lines) < 2:
            continue

        best_line = lines[1].strip().split()
        partition_number = best_line[0]

        partition_file = os.path.join(asap_dir, gene, f"{gene}.fasta.Partition_{partition_number}.csv")
        if not os.path.exists(partition_file):
            continue

        group_dict = {}
        with open(partition_file) as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) < 2:
                    continue
                strain, group = row[0].strip(), row[1].strip()
                group_dict[strain] = group

        matrix = build_partition_matrix(strains, group_dict)
        output_path = os.path.join(output_dir, f"{gene}.csv")
        save_partition_matrix(strains, matrix, output_path)
        print(f"Matrix saved for {gene}")

    except Exception as e:
        print(f"Error processing {gene}: {e}")

print("\n Partition matrices saved in:", output_dir)
