import os
import subprocess

input_dir = os.path.expanduser("~/core_genes_fasta")
output_dir = os.path.expanduser("~/core_genes_aligned")

os.makedirs(output_dir, exist_ok=True)

for fasta_file in os.listdir(input_dir):
    if fasta_file.endswith(".fasta"):
        input_path = os.path.join(input_dir, fasta_file)
        output_path = os.path.join(output_dir, fasta_file.replace(".fasta", "_aligned.fasta"))

        print(f"Aligning {fasta_file}...")
        subprocess.run(["mafft", "--auto", input_path], stdout=open(output_path, "w"))

print("All alignments saved to:", output_dir)
