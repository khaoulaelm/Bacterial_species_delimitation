#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")

import pandas as pd
import matplotlib.pyplot as plt
from statistics import multimode


def plateaus_by_value(df, col):
    thr = df["Threshold_int"].to_list()
    vals = df[col].to_list()

    ranges = []
    start = prev_thr = prev_val = None

    for t, v in zip(thr, vals):
        if prev_val is None:
            start, prev_thr, prev_val = t, t, v
            continue

        if (v == prev_val) and (t == prev_thr + 1):
            prev_thr = t
        else:
            ranges.append((prev_val, start, prev_thr))
            start, prev_thr, prev_val = t, t, v

    if prev_val is not None:
        ranges.append((prev_val, start, prev_thr))

    grouped = {}
    for v, a, b in ranges:
        grouped.setdefault(v, []).append((a, b))

    return grouped


def best_plateau(df, group_col):
    mode_g = multimode(df[group_col])[0]
    P = plateaus_by_value(df, group_col)
    ranges = P.get(mode_g, [])

    if not ranges:
        return mode_g, (None, None)

    # longest plateau
    best = sorted(ranges, key=lambda x: -(x[1] - x[0]))[0]
    return mode_g, best


def prepare_df(path):
    df = pd.read_csv(path)
    group_col = "Num_Groups" if "Num_Groups" in df.columns else "Num_VUB_Groups"

    df["Threshold"] = pd.to_numeric(df["Threshold"], errors="coerce")
    df = df.dropna(subset=["Threshold"]).copy()

    df["Threshold_int"] = df["Threshold"].round().astype(int)
    df = df.sort_values("Threshold_int").reset_index(drop=True)
    return df, group_col


def run(abgd_csv, asap_csv, selected_groups=11, outfile="ABGD_ASAP_thresholds_plateau.pdf"):
    abgd, abgd_col = prepare_df(abgd_csv)
    asap, asap_col = prepare_df(asap_csv)

    abgd_mode, abgd_plateau = best_plateau(abgd, abgd_col)
    asap_mode, asap_plateau = best_plateau(asap, asap_col)

    fig, ax = plt.subplots()

    # Curves
    ax.plot(abgd["Threshold"], abgd[abgd_col], color="blue", linewidth=2, label="ABGD VUB")
    ax.plot(asap["Threshold"], asap[asap_col], color="orange", linewidth=2, label="ASAP VUB")

    # Plateaus 
    a1, b1 = abgd_plateau
    if a1 is not None:
        ax.axvspan(a1, b1, alpha=0.15, color="blue", label="ABGD plateau")

    a2, b2 = asap_plateau
    if a2 is not None:
        ax.axvspan(a2, b2, alpha=0.15, color="orange", label="ASAP plateau")

    
    ax.axhline(selected_groups, linestyle="--", color="orange", linewidth=1)
    xmin = min(abgd["Threshold"].min(), asap["Threshold"].min())
    ax.text(xmin, selected_groups + 0.2, f"{selected_groups} groups", color="gray")

    ax.set_xlabel("Thresholds (conspecificity score)")
    ax.set_ylabel("Number of Groups")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(outfile, format="pdf", dpi=300)
    plt.close(fig)

    print("Saved:", outfile)
    print("ABGD best group (mode):", abgd_mode, " plateau:", abgd_plateau)
    print("ASAP best group (mode):", asap_mode, " plateau:", asap_plateau)


if __name__ == "__main__": 
    import sys 
    if len(sys.argv) == 3:
        abgd_csv = sys.argv[1] 
        asap_csv = sys.argv[2] 
    else: 
        abgd_csv = "/vub_abgdthreshold_summary.csv" 
        asap_csv = "/vub_asapthreshold_summary.csv" 
    run(abgd_csv, asap_csv)
