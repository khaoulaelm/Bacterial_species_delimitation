#!/usr/bin/env python3
"""
Author: Khaoula El Mchachti
Description: Scan conspecificity thresholds to evaluate how the number of inferred species groups changes with the number of core genes supporting strain grouping for ABGD and ASAP. For each threshold, a graph is constructed where strains are connected if their conspecificity score (number of genes assigning them to the same group) is greater than or equal to the threshold. The number of connected components corresponds to the number of inferred species groups. Plateaus (continuous ranges of thresholds producing the same number of groups) are detected, and the longest plateau is highlighted as the most stable species delimitation. Both methods are displayed on the same plot for easy comparison.
Input: ABGD_threshold_summary.csv and ASAP_threshold_summary.csv
Output: Combined PDF plot showing the number of groups vs thresholds for ABGD and ASAP
Date: 2026-03-30
Last modified: 2026-08-17
"""

import matplotlib
matplotlib.use("Agg")
import os
import pandas as pd
import matplotlib.pyplot as plt

# Find the directory containing this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Project root directory
PROJECT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "../.."))

def find_plateaus(df, ycol="Num_Groups"):

    thr = df["Threshold_int"].to_list()
    vals = df[ycol].to_list()

    plateaus = []
    
    start = prev_thr = prev_val = None

    for t, v in zip(thr, vals):
    
        if prev_val is None:
            start, prev_thr, prev_val = t, t, v
            continue

        if (v == prev_val) and (t == prev_thr + 1):
            prev_thr = t
            
        else:
            plateaus.append((prev_val, start, prev_thr))
            start, prev_thr, prev_val = t, t, v

    if prev_val is not None:
    
        plateaus.append((prev_val, start, prev_thr))
    
    return plateaus
    

def best_plateau(df, ycol="Num_Groups"):
    
    plateaus = find_plateaus(df, ycol)

    if not plateaus:
    
        return None, (None, None)

    # Select 1.longest plateau 2. Highest starting threshold in case of a tie
    def plateau_key(p):
        value, start, end = p
        length = end - start
        
        return(length, start)
        
    best = sorted(plateaus, key=plateau_key, reverse=True)[0]
    
    best_value, start, end = best
    
    return ( best_value, (start, end))


def prepare_df(path):

    df = pd.read_csv(path)

    df["Threshold"] = pd.to_numeric(df["Threshold"], errors="coerce")
    df = df.dropna(subset=["Threshold"]).copy()

    df["Threshold_int"] = df["Threshold"].round().astype(int)
    df = df.sort_values("Threshold_int").reset_index(drop=True)
    
    return df

def plot_combined( abgd, asap, abgd_val, abgd_range, asap_val, asap_range, outfile):
    fig, ax = plt.subplots()
    
    # ABGD curve
    ax.plot(abgd["Threshold_int"], abgd["Num_Groups"], color="blue", linewidth=2, label="ABGD VUB")
    # ABGD best plateau
    a1, b1 = abgd_range

    if a1 is not None:
        ax.axvspan(a1, b1, alpha=0.15, color="blue", label=f"ABGD plateau {a1}-{b1}")
    
    # ASAP curve
    ax.plot(asap["Threshold_int"], asap["Num_Groups"], color="orange", linewidth=2, label="ASAP VUB")
    # ASAP best plateau
    a2, b2 = asap_range

    if a2 is not None:
        ax.axvspan(a2, b2, alpha=0.15, color="orange", label=f"ASAP plateau {a2}-{b2}")  
        
    # Figure labels
    ax.set_xlabel("Thresholds (conspecificity score)")
    ax.set_ylabel("Number of Groups")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(outfile, format="pdf", dpi=300)
    plt.close(fig) 
    
def run():
    
    # ABGD input
    abgd_csv = os.path.join(
        PROJECT_DIR,
        "3_species_delimitation_methods",
        "4_CGCD_approach",
        "ABGD_threshold_scan",
        "ABGD_threshold_summary.csv"
    )
    
    # ASAP input
    asap_csv = os.path.join(
        PROJECT_DIR,
        "3_species_delimitation_methods",
        "4_CGCD_approach",
        "ASAP_threshold_scan",
        "ASAP_threshold_summary.csv"
    )    
    
    # Read both datasets
    abgd = prepare_df(abgd_csv)
    asap = prepare_df(asap_csv)
    
    # Detect best plateau
    abgd_val, abgd_range = best_plateau(abgd, ycol="Num_Groups")
    asap_val, asap_range = best_plateau(asap, ycol="Num_Groups")
    
    # Output directories

    combined_out_dir = os.path.join(
        PROJECT_DIR,
        "3_species_delimitation_methods",
        "4_CGCD_approach",
        "combined_plots"
    )

    os.makedirs(combined_out_dir, exist_ok=True)    
   
    # Combined plot

    combined_pdf = os.path.join(
        combined_out_dir,
        "Groups_vs_threshold_abgd_asap.pdf"
    )

    plot_combined(
        abgd,
        asap,
        abgd_val,
        abgd_range,
        asap_val,
        asap_range,
        combined_pdf
    )

    print( "Combined plot saved to:", combined_pdf)

if __name__ == "__main__":

    run()    
