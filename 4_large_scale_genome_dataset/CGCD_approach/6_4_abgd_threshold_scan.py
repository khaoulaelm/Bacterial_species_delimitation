#!/usr/bin/env python3
"""
Author: Khaoula El Mchachti
Description: This script analyzes an ABGD conspecificity matrix to explore how strain groupings change across thresholds, where the threshold represents the number of shared core genes supporting the grouping. In this analysis, only high thresholds (>= 80% of the maximum value) are considered. For each threshold, a network is constructed where strains are connected if their conspecificity score meets or exceeds the threshold, and groups are defined as connected components.
The script calculates both the total number of groups (all strains) and the number of groups 
containing at least one VUB strain. The results are saved in a summary file,  allowing visualization of how grouping patterns vary across the selected threshold range.
Input: ABGD_conspecificity_matrix.csv, VUBstrains.csv
Output: ABGD_groupings_all
Date: 2026-04-14
""" 

import pandas as pd
import networkx as nx
import os

# Input files
matrix_path = "ABGD_conspecificity_matrix/ABGD_conspecificity_matrix.csv"
vub_strains_path = "VUBstrains.csv"

# Load data
df = pd.read_csv(os.path.expanduser(matrix_path), index_col=0)
vub_df = pd.read_csv(os.path.expanduser(vub_strains_path))
vub_strains = set(vub_df["Strain"].tolist())  

# Output directory 
output_dir = "ABGD_groupings_all"
os.makedirs(output_dir, exist_ok=True)
summary_path = os.path.join(output_dir, "ABGD_groupings_summary.csv")

# Scan thresholds
vals = df.values[df.values > 0]
min_threshold = int(vals.min())
max_threshold = int(vals.max())

# Set 80% cutoff
threshold_80 = int(0.8 * max_threshold)

print(f"Scanning thresholds from {threshold_80} to {max_threshold}")

strains = list(df.index)
n = len(strains)
group_summary = []

for threshold in range(threshold_80, max_threshold + 1):

    G = nx.Graph()
    G.add_nodes_from(df.index)

    # Add edges where value >= threshold
    for i, strain_i in enumerate(df.index):
        for j, strain_j in enumerate(df.columns):
            if i < j and df.loc[strain_i, strain_j] >= threshold:
                G.add_edge(strain_i, strain_j)

    # Connected components = groups
    components = list(nx.connected_components(G))
    total_num_groups = len(components)

    # Count how many components contain at least one VUB strain
    num_vub_groups = sum(1 for comp in components if any(s in vub_strains for s in comp))

    group_summary.append((threshold, num_vub_groups, total_num_groups))

# Save summary 
summary_df = pd.DataFrame(
    group_summary,
    columns=["Threshold", "Num_VUB_Groups", "Num_Total_Groups"]
)
summary_path = os.path.join(output_dir, "ABGD_groupings_summary.csv")
summary_df.to_csv(summary_path, index=False)

print(f"Summary saved to: {summary_path}")
