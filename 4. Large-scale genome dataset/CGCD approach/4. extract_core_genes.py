import os
import pandas as pd
from Bio import SeqIO

# Paths
csv_path = os.path.expanduser("~/roary_results/gene_presence_absence.csv")
prokka_base_dir = os.path.expanduser("~/prokka_results")
output_dir = os.path.expanduser("~/core_genes_fasta3")
missing_log_path = os.path.join(output_dir, "missing_genes_log.csv")

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load Roary gene matrix 
df = pd.read_csv(csv_path, low_memory=False)
strain_columns = df.columns[14:]

# Select strict core genes (present in ALL strains)
core_genes = df[df["No. isolates"].astype(float) == len(strain_columns)]
print(f" Found {len(core_genes)} strict core genes (present in all {len(strain_columns)} strains).")

# Initialize missing data logger 
missing_data = []

# Loop over core genes 
for idx, row in core_genes.iterrows():
    gene_name = row["Gene"] if pd.notna(row["Gene"]) else f"group_{idx+1}"
    output_fasta = os.path.join(output_dir, f"{gene_name}.fasta")
    sequences_written = 0

    with open(output_fasta, "w") as out_fasta:
        for strain in strain_columns:
            gene_ids_str = row[strain]

            if pd.isna(gene_ids_str):
                missing_data.append((gene_name, strain, "Missing gene ID in Roary"))
                continue

            # Split IDs by semicolon or spaces
            gene_ids = [g.strip() for g in str(gene_ids_str).replace(";", " ").split() if g.strip()]

            strain_ffn = os.path.join(prokka_base_dir, strain, f"{strain}.ffn")
            if not os.path.isfile(strain_ffn):
                missing_data.append((gene_name, strain, "Missing .ffn file"))
                continue

            # Load ALL records (list to keep duplicates)
            records = list(SeqIO.parse(strain_ffn, "fasta"))
            matched = False

            for gene_id in gene_ids:
                for record in records:
                    if gene_id in record.id:  # Match ONLY by ID
                        record.id = strain
                        record.description = ""
                        SeqIO.write(record, out_fasta, "fasta")
                        print(f"{gene_name}  {strain} (selected {record.id})")
                        sequences_written += 1
                        matched = True
                        break
                if matched:
                    break

            if not matched:
                missing_data.append((gene_name, strain, "Gene ID not found in .ffn"))

    # Delete FASTA if no sequences written for this gene
    if sequences_written == 0:
        os.remove(output_fasta)

# Save missing gene log 
pd.DataFrame(missing_data, columns=["Gene", "Strain", "Reason"]).to_csv(missing_log_path, index=False)
print(f"Missing gene log saved to: {missing_log_path}")

print(f"\n Extraction complete. FASTA files saved in: {output_dir}")
