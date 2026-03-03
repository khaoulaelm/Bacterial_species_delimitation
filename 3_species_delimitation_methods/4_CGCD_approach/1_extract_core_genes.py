import os
import pandas as pd
from Bio import SeqIO

# Paths
csv_path = os.path.expanduser("~/roary/gene_presence_absence.csv")
prokka_base_dir = os.path.expanduser("~/prokka")
output_dir = os.path.expanduser("~/core_genes_fasta")

# Make sure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load CSV
df = pd.read_csv(csv_path)

# Get strain names from column headers
strain_columns = df.columns[14:]  

# Filter core genes (present in all 47 isolates)
core_genes = df[df["No. isolates"] == 47]

print(f"Found {len(core_genes)} core genes")

# Iterate through each core gene
for idx, row in core_genes.iterrows():
    gene_name = row["Gene"] if pd.notna(row["Gene"]) else f"group_{idx+1}"
    output_fasta = os.path.join(output_dir, f"{gene_name}.fasta")

    with open(output_fasta, "w") as out_fasta:
        for strain in strain_columns:
            gene_ids = row[strain]
            if pd.isna(gene_ids):
                continue

            for gene_id in str(gene_ids).split(';'):
                gene_id = gene_id.strip()
                if not gene_id:
                    continue

                strain_ffn_path = os.path.join(prokka_base_dir, strain, f"{strain}.ffn")
                if not os.path.isfile(strain_ffn_path):
                    print(f" Missing file for {strain}: {strain_ffn_path}")
                    continue

                for record in SeqIO.parse(strain_ffn_path, "fasta"):
                    if record.id == gene_id:
                        record.id = strain  # Rename to strain name
                        record.description = ""
                        SeqIO.write(record, out_fasta, "fasta")
                        break

print("Extraction complete. FASTA files saved in:", output_dir)
