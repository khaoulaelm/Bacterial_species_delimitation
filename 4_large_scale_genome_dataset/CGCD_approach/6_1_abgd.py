#!/usr/bin/env python3

"""
Author: Khaoula El Mchachti
Description: Run ABGD on each aligned core-gene (one gene at a time)
Input: core_genes_aligned/
Output: ABGD_results/<gene_name>/ (ABGD results per gene), failed_abgd_genes.txt (log file listing failures) and abgd_warnings.txt (listing genes with warnings, non-zero exit codes, or stderr output)
Date: 2026-04-12
"""

import os, glob, subprocess

abgd_exec = os.path.expanduser("~/Bacterial_species_delimitation/3_species_delimitation_methods/ABGD/abgd")
input_dir = os.path.expanduser("core_genes_aligned")
main_output_dir = os.path.expanduser("ABGD_results")
failed_genes_file = os.path.expanduser("ABGD_results/failed_abgd_genes.txt")
warnings_file = os.path.join(main_output_dir, "abgd_warnings.txt")

os.makedirs(main_output_dir, exist_ok=True)

def abgd_outputs_exist(outdir, basename):
    # any partition file or the result cvs counts as success
    parts = glob.glob(os.path.join(outdir, f"{basename}.part*.txt"))
    cvs   = os.path.join(outdir, f"{basename}.res.cvs")
    return (len(parts) > 0) or os.path.isfile(cvs)

with open(failed_genes_file, "w") as failed, open(warnings_file, "w") as warnlog:
    failed.write("Failed genes:\n")
    for fname in os.listdir(input_dir):
        if not fname.endswith(".fasta"):
            continue

        gene_path = os.path.join(input_dir, fname)
        gene_name = os.path.splitext(fname)[0]
        gene_output_dir = os.path.join(main_output_dir, gene_name)
        os.makedirs(gene_output_dir, exist_ok=True)

        cmd = [abgd_exec, "-a", "-d", "JC69" , "-o", gene_output_dir, gene_path]
        print(f"Running ABGD on: {gene_name}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Save raw logs for debugging
        with open(os.path.join(gene_output_dir, "abgd.stdout.txt"), "w") as f: f.write(result.stdout or "")
        with open(os.path.join(gene_output_dir, "abgd.stderr.txt"), "w") as f: f.write(result.stderr or "")

        outputs_ok = abgd_outputs_exist(gene_output_dir, gene_name)

        if outputs_ok:
            # success, but keep a heads-up if something looked noisy
            if result.returncode != 0 or ("error" in (result.stderr or "").lower()):
                warnlog.write(f"[{gene_name}] returncode={result.returncode}\n")
                if result.stderr:
                    warnlog.write(result.stderr + "\n---\n")
            # Optional: tag “single-partition” if CVS is empty
            cvs_path = os.path.join(gene_output_dir, f"{gene_name}.res.cvs")
            if os.path.isfile(cvs_path):
                try:
                    # if only header or blank, likely single partition
                    lines = [ln for ln in open(cvs_path) if ln.strip()]
                    if len(lines) <= 1:
                        print(f"ABGD completed (single partition) for {gene_name}.")
                    else:
                        print(f"ABGD completed successfully for {gene_name}.")
                except Exception:
                    print(f"ABGD completed for {gene_name} (could not read CVS).")
            else:
                print(f"ABGD completed successfully for {gene_name}.")
        else:
            print(f"ABGD produced no outputs for {gene_name}. Marking as failed.")
            failed.write(f"{gene_name}\n")

print("All genes processed with ABGD.")
print("Failed genes are logged in:", failed_genes_file)
print("Any warnings (non-zero exit or stderr noise) are in:", warnings_file)
print("\nNote: Check failed and warning logs for genes affected by ABGD issues (often memory-related).")
