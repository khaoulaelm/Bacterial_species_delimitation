#!/usr/bin/env python3
"""
Author: Khaoula El Mchachti
Description: Scan conspecificity thresholds to evaluate how the number of inferred species groups changes with the number of core genes supporting strain grouping. For each threshold, a graph is constructed where strains are connected if their conspecificity score (number of genes assigning them to the same group) is greater than or equal to the threshold. The number of connected components corresponds to the number of inferred species groups. This analysis focuses specifically on VUB strain groups rather than all strains. Plateaus (continuous ranges of thresholds producing the same number of groups) are detected, and the longest plateau is highlighted as the most stable species delimitation.
Input: ABGD_groupings_summary.csv
Output:Groups_vs_thresholds_abgd.pdf
Date: 2026-04-14
"""

import matplotlib
matplotlib.use("Agg")
import os
import pandas as pd
import matplotlib.pyplot as plt
from statistics import multimode

def find_plateaus(df, ycol="Num_VUB_Groups"):
    
    thr = df["Threshold_int"].to_list()
    vals = df[ycol].to_list()

    plateaus = []
    start = prev_thr = prev_val = None

    for t, v in zip(thr, vals):
        if prev_val is None:
            start, prev_thr, prev_val = t, t, v
            continue

        # Continue plateau if same value and thresholds are consecutive
        if v == prev_val and t == prev_thr + 1:
            prev_thr = t
        else:
            # Close previous plateau
            plateaus.append((prev_val, start, prev_thr))
            # Start new plateau
            start, prev_thr, prev_val = t, t, v

    # Close last plateau
    if prev_val is not None:
        plateaus.append((prev_val, start, prev_thr))

    return plateaus

def best_plateau(df, ycol="Num_VUB_Groups"):

    plateaus = find_plateaus(df, ycol=ycol)
    if not plateaus:
        return None, (None, None), []

    # Sort by: longest length first, then highest start threshold (higher gene agreement)
    
    def plateau_key(p):
        val, a, b = p
        length = b - a
        return (length, a)

    best = sorted(plateaus, key=plateau_key, reverse=True)[0]
    best_val, a, b = best
    return best_val, (a, b), plateaus
    
def prepare_df(path):
    df = pd.read_csv(path)

    # Require expected columns
    if "Threshold" not in df.columns or "Num_VUB_Groups" not in df.columns:
        raise ValueError("CSV must contain columns: 'Threshold' and 'Num_VUB_Groups'.")

    df["Threshold"] = pd.to_numeric(df["Threshold"], errors="coerce")
    df = df.dropna(subset=["Threshold"]).copy()

    # Create integer thresholds for plateau detection
    df["Threshold_int"] = df["Threshold"].round().astype(int)

    df = df.sort_values("Threshold_int").reset_index(drop=True)
    return df

def plot_and_save(df, ycol, best_val, best_range, outfile, line_color="blue"):
    fig, ax = plt.subplots()

    ax.plot(df["Threshold_int"], df[ycol], linewidth=2, color=line_color)

    a, b = best_range
    if best_val is not None:
        ax.axhline(best_val, linestyle="--", color=line_color,
                   label=f"Best plateau: {best_val} groups")

    if a is not None:
        ax.axvspan(a, b, alpha=0.15, color=line_color, label=f"Plateau {a}-{b}")

    ax.set_xlabel("Thresholds (conspecificity score)")
    ax.set_ylabel("Number of groups")
    ax.legend()
    fig.tight_layout()

    fig.savefig(outfile, format="pdf", dpi=300)
    plt.close(fig)


def run(csv_path, prefix="Groups_vs_threshold"):
    df = prepare_df(csv_path)

    best_val, best_range, plateaus = best_plateau(df, ycol="Num_VUB_Groups")

    out_dir= "ABGD_plots"
    os.makedirs(out_dir, exist_ok=True)

    out_pdf = os.path.join(out_dir, f"{prefix}.pdf")
    plot_and_save(df, "Num_VUB_Groups", best_val, best_range, out_pdf, line_color="blue")

    print("Saved:", out_pdf)
    if best_val is not None:
        a, b = best_range
        print(f"Best plateau = {best_val} groups (Threshold {a}–{b})")
    else:
        print("No plateau detected.")


if __name__ == "__main__":
    import sys

    if len(sys.argv) >= 2:
        csv_path = sys.argv[1]
    else:
        csv_path = "ABGD_groupings_all/ABGD_groupings_summary.csv"

    run(csv_path, prefix="Groups_vs_thresholds_abgd")
