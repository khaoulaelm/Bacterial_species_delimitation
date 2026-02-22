import os
import glob
import subprocess

# Paths
asap_exec   = os.path.expanduser("~/ASAP/asap")
input_dir   = os.path.expanduser("~/core_genes_aligned")
output_dir  = os.path.expanduser("~/ASAP_results")
failed_log  = os.path.join(output_dir, "failed_asap_genes.txt")
warnings_log = os.path.join(output_dir, "asap_warnings.txt")

os.makedirs(output_dir, exist_ok=True)

def find_spart(out_dir, gene_base):
    """
    Return the best-guess .spart path if present, else None.
    We try the most common names:
      gene_base.spart
      gene_base*.spart
    """
    candidates = [
        os.path.join(out_dir, f"{gene_base}.spart"), 
    ] + glob.glob(os.path.join(out_dir, f"{gene_base}*.spart"))
    for p in candidates:
        if os.path.isfile(p):
            return p
    return None

def spart_has_partitions(spart_path):
    """Quick check: does the .spart file contain any partition lines?"""
    try:
        with open(spart_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                # ASAP .spart has blocks like: "Partition X;" followed by assignments
                if line.strip().lower().startswith("partition"):
                    return True
        # Fallback: non-empty file counts as something
        return os.path.getsize(spart_path) > 0
    except Exception:
        return False

# Initialize logs
with open(failed_log, "w") as flog:
    flog.write("Failed Genes:\n")
with open(warnings_log, "w") as wlog:
    wlog.write("Warnings (non-zero return code or stderr present):\n\n")

# Loop through aligned files
for fname in os.listdir(input_dir):
    if not fname.endswith(".fasta"):
        continue

    gene_base = os.path.splitext(fname)[0]          
    in_path   = os.path.join(input_dir, fname)
    out_dir   = os.path.join(output_dir, gene_base)
    os.makedirs(out_dir, exist_ok=True)

    # Skip if already processed (a .spart exists)
    existing_spart = find_spart(out_dir, gene_base)
    if existing_spart:
        print(f"Skipping {gene_base} (already processed: {os.path.basename(existing_spart)})")
        continue

    print(f"\nRunning ASAP for {gene_base}...")
    try:
        result = subprocess.run(
            [asap_exec, "-a", "-o", out_dir, in_path],
            capture_output=True,
            text=True,
            cwd=out_dir
            # , timeout=90
        )

        # Save raw logs for debugging
        with open(os.path.join(out_dir, "asap.stdout.txt"), "w") as f: f.write(result.stdout or "")
        with open(os.path.join(out_dir, "asap.stderr.txt"), "w") as f: f.write(result.stderr or "")

        # Decide success based on outputs
        spart_path = find_spart(out_dir, gene_base)
        success = bool(spart_path and spart_has_partitions(spart_path))

        if success:
            # It succeeded; still record warnings if something looked off
            if result.returncode != 0 or (result.stderr and result.stderr.strip()):
                with open(warnings_log, "a") as wlog:
                    wlog.write(f"[{gene_base}] returncode={result.returncode}\n")
                    if result.stderr:
                        wlog.write(result.stderr + "\n---\n")
            # Optional: detect “single partition” by counting blocks
            n_parts = 0
            if spart_path:
                with open(spart_path, "r", encoding="utf-8", errors="ignore") as f:
                    n_parts = sum(1 for ln in f if ln.strip().lower().startswith("partition"))
            if n_parts <= 1:
                print(f"ASAP completed (single partition) for {gene_base}.")
            else:
                print(f"ASAP completed for {gene_base} ({n_parts} partitions).")
        else:
            print(f"ASAP produced no usable .spart for {gene_base}. Marking as failed.")
            with open(failed_log, "a") as flog:
                flog.write(f"{gene_base}\n")

    except subprocess.TimeoutExpired:
        print(f"Timeout: ASAP took too long on {gene_base}. Skipping.")
        with open(failed_log, "a") as flog:
            flog.write(f"{gene_base} (timeout)\n")

    except Exception as e:
        print(f"Exception for {gene_base}: {e}")
        with open(failed_log, "a") as flog:
            flog.write(f"{gene_base} (exception)\n")

print("\nAll genes processed.")
print("Failed genes listed in:", failed_log)
print("Warnings (non-zero exit / stderr) in:", warnings_log)
