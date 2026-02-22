#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")  
import pandas as pd, math
import matplotlib.pyplot as plt
from statistics import multimode, median
from matplotlib.backends.backend_pdf import PdfPages


def plateaus_by_value(df, col="Num_VUB_Groups"):
    ranges, start, prev_thr, prev_val = [], None, None, None
    for thr, val in zip(df["Threshold"], df[col]):
        if prev_val is None:
            start, prev_thr, prev_val = thr, thr, val
            continue
        if (val == prev_val) and (thr == prev_thr + 1):
            prev_thr = thr
        else:
            ranges.append((prev_val, start, prev_thr))
            start, prev_thr, prev_val = thr, thr, val
    ranges.append((prev_val, start, prev_thr))
    grouped = {}
    for v, a, b in ranges:
        grouped.setdefault(v, []).append((a, b))
    return grouped

def plateau_mode_and_prop(df, a, b):
    sub = df[(df["Threshold"] >= a) & (df["Threshold"] <= b)]
    m = multimode(sub["Num_Total_Groups"])[0]
    prop = (sub["Num_Total_Groups"] == m).mean()
    return m, float(prop)

def plateau_length(ranges):
    return int(sum(b - a + 1 for a, b in ranges))

def choose_best_plateau_for_g(df, ranges_g):
    scored = []
    for a, b in ranges_g:
        m, p = plateau_mode_and_prop(df, a, b)
        scored.append((p, b - a + 1, a, b, m))
    scored.sort(reverse=True)
    p, length, a, b, m = scored[0]
    return (a, b), m, p

def best_group_number(abgd_df, asap_df, kappa=20):
    P_abgd = plateaus_by_value(abgd_df, "Num_VUB_Groups")
    P_asap = plateaus_by_value(asap_df, "Num_VUB_Groups")
    candidates = sorted(set(P_abgd.keys()) & set(P_asap.keys()))
    best = None
    for g in candidates:
        La = plateau_length(P_abgd[g])
        Ls = plateau_length(P_asap[g])
        (aA, bA), Ta, pa = choose_best_plateau_for_g(abgd_df, P_abgd[g])
        (aS, bS), Ts, ps = choose_best_plateau_for_g(asap_df, P_asap[g])
        score = (La * Ls) * (pa * ps) * math.exp(-abs(Ta - Ts) / 20)
        info = {
            "g": int(g),
            "score": float(score),
            "ABGD": {"plateau": (int(aA), int(bA)), "total_mode": int(Ta), "prop": float(pa), "sum_length": int(La)},
            "ASAP": {"plateau": (int(aS), int(bS)), "total_mode": int(Ts), "prop": float(ps), "sum_length": int(Ls)},
        }
        if (best is None) or (score > best["score"]):
            best = info
    return best

def pick_threshold_in_plateau(df, plateau, target_g):
    a, b = plateau
    sub = df[(df["Threshold"] >= a) & (df["Threshold"] <= b)]
    T_mode = multimode(sub["Num_Total_Groups"])[0]
    cand = sub[(sub["Num_VUB_Groups"] == target_g) & (sub["Num_Total_Groups"] == T_mode)]
    if not cand.empty:
        thr = int(round(median(cand["Threshold"].astype(int))))
    else:
        thr = int(round((a + b) / 2))
    if (df["Threshold"] == thr).any():
        row = df[df["Threshold"] == thr].iloc[0]
    else:
        idx = (df["Threshold"] - thr).abs().idxmin()
        row = df.loc[idx]
        thr = int(row["Threshold"])
    return {
        "threshold": thr,
        "vub_groups": int(row["Num_VUB_Groups"]),
        "total_groups": int(row["Num_Total_Groups"]),
        "total_mode": int(T_mode),
    }

def run(abgd_csv, asap_csv, out_prefix="best_groups"):
    abgd = pd.read_csv(abgd_csv)
    asap = pd.read_csv(asap_csv)

    best = best_group_number(abgd, asap, kappa=20)
    g_star = best["g"]

    abgd_pick = pick_threshold_in_plateau(abgd, best["ABGD"]["plateau"], g_star)
    asap_pick = pick_threshold_in_plateau(asap, best["ASAP"]["plateau"], g_star)

    # COMBINED PLOT 
    from matplotlib.backends.backend_pdf import PdfPages
    pdf_file = f"{out_prefix}_combined.pdf"
    with PdfPages(pdf_file) as pdf:

        plt.figure(figsize=(8,5))

        # ABGD curve
        plt.plot(abgd["Threshold"], abgd["Num_VUB_Groups"], label="ABGD VUB", color="blue")
        aA, bA = best["ABGD"]["plateau"]
        plt.axvspan(aA, bA, alpha=0.1, color="blue", label=f"ABGD plateau")
        plt.axhline(abgd_pick["vub_groups"], linestyle="--", color="blue")
        # annotate ABGD VUB
        plt.text(abgd["Threshold"].min(), abgd_pick["vub_groups"] + 0.3,  # slightly above line
        f"{abgd_pick['vub_groups']} groups", color="gray", fontsize=9)

        

        # ASAP curve
        plt.plot(asap["Threshold"], asap["Num_VUB_Groups"], label="ASAP VUB", color="orange")
        aS, bS = best["ASAP"]["plateau"]
        plt.axvspan(aS, bS, alpha=0.1, color="orange", label=f"ASAP plateau")
        plt.axhline(asap_pick["vub_groups"], linestyle="--", color="orange")

        plt.xlabel("Thresholds (number of core-genes)")
        plt.ylabel("Number of Groups")

        plt.legend()
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    print(f"Combined PDF saved to {pdf_file}")

    # CSV SUMMARY 
    summary = pd.DataFrame([
        {"method":"ABGD","chosen_g":g_star,"threshold":abgd_pick["threshold"],
         "vub_groups":abgd_pick["vub_groups"],"total_groups":abgd_pick["total_groups"],
         "plateau_start":best["ABGD"]["plateau"][0],"plateau_end":best["ABGD"]["plateau"][1],
         "total_mode_used":best["ABGD"]["total_mode"],"mode_proportion":best["ABGD"]["prop"],
         "sum_plateau_length":best["ABGD"]["sum_length"]},
        {"method":"ASAP","chosen_g":g_star,"threshold":asap_pick["threshold"],
         "vub_groups":asap_pick["vub_groups"],"total_groups":asap_pick["total_groups"],
         "plateau_start":best["ASAP"]["plateau"][0],"plateau_end":best["ASAP"]["plateau"][1],
         "total_mode_used":best["ASAP"]["total_mode"],"mode_proportion":best["ASAP"]["prop"],
         "sum_plateau_length":best["ASAP"]["sum_length"]}
    ])
    summary.to_csv(f"{out_prefix}_summary.csv", index=False)
    print(summary.to_string(index=False))
    print("Chosen g:", g_star)
    return g_star


if __name__ == "__main__":
    # Example usage:
    # python plot_best_groups.py vub_ABGDthreshold_summary.csv vub_ASAPthreshold_summary.csv
    import sys
    abgd_csv = sys.argv[1] if len(sys.argv) > 1 else "/ABGD_groupings/ABGD_groupings_summary.csv"
    asap_csv = sys.argv[2] if len(sys.argv) > 2 else "/ASAP_groupings/ASAP_groupings_summary.csv"
    run(abgd_csv, asap_csv)
