#!/usr/bin/env python3
"""
Author: Khaoula El Mchachti
Description: Extract strain groups for each threshold within a detected plateau of the conspecificity threshold scan.
Input: conspecificity_matrix.csv
Output: ABGD_groups_plateau/groups_t{threshold}.csv
Date: 2026-03-20
"""
import pandas as pd
import networkx as nx
import os

# Load conspecificity matrix
df = pd.read_csv("ABGD_conspecificity_matrix/ABGD_conspecificity_matrix.csv", index_col=0)

# All strains
strains = df.index.tolist()

# Plateau range (from your plot)
print("NOTE: Please update START and END to match the plateau range from your plot before running.")
START = 1494
END = 2026

# Output directory
output_dir = "ABGD_groups_plateau"
os.makedirs(output_dir, exist_ok=True)

print("Extracting groups for thresholds:", START, "to", END)


for threshold in range(START, END + 1):

    G = nx.Graph()
    G.add_nodes_from(strains)

    # Add edges if value >= threshold
    for i, s1 in enumerate(strains):
        for j, s2 in enumerate(strains):
            if i < j and df.loc[s1, s2] >= threshold:
                G.add_edge(s1, s2)

    # Connected components = groups
    components = list(nx.connected_components(G))

    # Save groups
    rows = []
    for gid, members in enumerate(components, start=1):
        for strain in members:
            rows.append({
                "Strain": strain,
                "Group": gid
            })

    pd.DataFrame(rows).to_csv(
        f"ABGD_groups_plateau/groups_t{threshold}.csv",
        index=False
    )

print("Done. Groups extracted for all thresholds in the plateau.")
print(f"Output directory: {os.path.abspath(output_dir)}")
