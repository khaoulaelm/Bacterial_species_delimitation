#!/usr/bin/env python3
import os
import csv

# Base directory for ASAP results
ASAP_DIR = os.path.expanduser("~/ASAP_results")
OUT_DIR  = os.path.expanduser("~/asap_partition_matrices")
os.makedirs(OUT_DIR, exist_ok=True)

def read_selected_partition_number(res_cvs_path: str) -> str | None:
    """Return the selected partition number from *.res.cvs (always it's the first partition)."""
    try:
        with open(res_cvs_path, "r") as f:
            lines = [ln.strip() for ln in f.readlines()]
        if len(lines) < 2:
            return None
        return lines[1].split()[0]
    except Exception:
        return None

def load_partition_map(partition_csv_path: str) -> dict[str, str]:
    """
    Load strain -> group mapping from the ASAP partition CSV (two columns: strain, group).
    Returns a dict {strain: group}.
    """
    m = {}
    with open(partition_csv_path, "r") as f:
        rdr = csv.reader(f)
        for row in rdr:
            if len(row) < 2:
                continue
            strain = row[0].strip()
            group  = row[1].strip()
            if strain:
                m[strain] = group
    return m

def build_matrix(strains: list[str], group_map: dict[str, str]) -> list[list[int]]:
    """
    Create a symmetric 0/1 matrix over the provided strain list.
    """
    # Assign a "missing" group label for absent strains so diagonal stays 1.
    MISS = "__MISSING__"
    lbl = [group_map.get(s, MISS) for s in strains]
    n = len(strains)
    M = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            M[i][j] = int(lbl[i] == lbl[j])
    return M

def save_matrix(strains: list[str], mat: list[list[int]], out_csv: str):
    """Write header + square matrix with row/col labels."""
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + strains)
        for s, row in zip(strains, mat):
            w.writerow([s] + row)

def main():
    # each subfolder of ASAP_DIR is a gene folder
    genes = [d for d in os.listdir(ASAP_DIR) if os.path.isdir(os.path.join(ASAP_DIR, d))]
    if not genes:
        print("No gene directories found in:", ASAP_DIR)
        return

    for gene in sorted(genes):
        gene_dir = os.path.join(ASAP_DIR, gene)
        res_cvs  = os.path.join(gene_dir, f"{gene}.fasta.res.cvs")
        if not os.path.exists(res_cvs):
            print(f"[{gene}] Missing {gene}.fasta.res.cvs — skipping.")
            continue

        part_num = read_selected_partition_number(res_cvs)
        if not part_num:
            print(f"[{gene}] Could not read selected partition number — skipping.")
            continue

        part_csv = os.path.join(gene_dir, f"{gene}.fasta.Partition_{part_num}.csv")
        if not os.path.exists(part_csv):
            print(f"[{gene}] Partition file not found: {os.path.basename(part_csv)} — skipping.")
            continue

        group_map = load_partition_map(part_csv)
        if not group_map:
            print(f"[{gene}] Partition file empty or invalid — skipping.")
            continue

        strains = sorted(group_map.keys())  # use only strains present in this partition file
        mat = build_matrix(strains, group_map)
        out_csv = os.path.join(OUT_DIR, f"{gene}.csv")
        save_matrix(strains, mat, out_csv)
        print(f"[{gene}] Matrix saved {out_csv}")

    print("\nAll ASAP partition matrices generated. Output:", OUT_DIR)

if __name__ == "__main__":
    main()
