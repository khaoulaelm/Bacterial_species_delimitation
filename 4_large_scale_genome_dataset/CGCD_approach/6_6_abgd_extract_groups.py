#!/usr/bin/env python3
"""
Author: Khaoula El Mchachti
Description: Extract strain groups for each threshold within a detected plateau of the conspecificity threshold scan, for both all strains and the VUB strains.
Input: ABGD_conspecificity_matrix.csv
Output: ABGD_groups_extraction
Date: 2026-04-14
"""

import pandas as pd
import networkx as nx
import os


# Input files
matrix_csv = "ABGD_conspecificity_matrix/ABGD_conspecificity_matrix.csv"   
vub_csv    = "VUBstrains.csv"        


# Load data
df = pd.read_csv(matrix_csv, index_col=0)
strains = df.index.tolist()

vub_df = pd.read_csv(vub_csv)
vub_strains = set(vub_df["Strain"].tolist())


# Manual plateau you choose 
print("NOTE: Please update START and END to match the plateau range from your plot before running.")
START = 185       
END   = 214

# Output directory 
out_dir = "ABGD_groups_extraction"
os.makedirs(out_dir, exist_ok=True)

print("Extracting groups for thresholds:", START, "to", END)

# Loop through plateau 
for threshold in range(START, END + 1):

    G = nx.Graph()
    G.add_nodes_from(strains)

    # Add edges if value >= threshold
    for i in range(len(strains)):
        s1 = strains[i]
        for j in range(i + 1, len(strains)):
            s2 = strains[j]
            if df.loc[s1, s2] >= threshold:
                G.add_edge(s1, s2)

    # Connected components = ABGD groups
    components = list(nx.connected_components(G))

    # Save ALL strain assignments (group IDs based on full components)
    all_rows = []
    for gid, members in enumerate(components, start=1):
        for strain in members:
            all_rows.append({"Strain": strain, "Group": gid})

    all_df = pd.DataFrame(all_rows).sort_values(["Group", "Strain"])
    all_out = os.path.join(out_dir, f"ABGD_groups_t{threshold}_ALL.csv")
    all_df.to_csv(all_out, index=False)

    # Save VUB-only assignments (same group IDs as ALL)
    vub_rows = []
    for gid, members in enumerate(components, start=1):
        for strain in members:
            if strain in vub_strains:
                vub_rows.append({"Strain": strain, "Group": gid})

    vub_df_out = pd.DataFrame(vub_rows).sort_values(["Group", "Strain"])
    vub_out = os.path.join(out_dir, f"ABGD_groups_t{threshold}_VUB.csv")
    vub_df_out.to_csv(vub_out, index=False)

    # Quick log
    num_total_groups = len(components)
    num_vub_groups = vub_df_out["Group"].nunique() if not vub_df_out.empty else 0
    print(f"t={threshold}: total_groups={num_total_groups}, vub_groups_present={num_vub_groups}")

print("Done. Saved ALL + VUB group files for all thresholds in the plateau.")
