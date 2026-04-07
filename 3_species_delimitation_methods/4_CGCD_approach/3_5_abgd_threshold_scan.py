#!/usr/bin/env python3
"""
Author: Khaoula El Mchachti
Description: Scan conspecificity thresholds to determine how the number of species groups changes based on the number of shared core genes supporting the grouping. The threshold corresponds to the minimum number of core genes that must assign two strains to the same group. For each threshold, a graph is built where strains are connected if their conspecificity score (number of genes supporting their grouping) is ≥ threshold. The number of connected components represents the number of inferred species groups.
Remark: The CGCD approach is fundamentally threshold-free, as species boundaries can be inferred by examining how the number of groups changes across the entire range of thresholds. However, for visualization purposes, it is often useful to focus on the region where a high proportion of genes agree on the grouping (e.g., >50% of core genes). Users may first inspect the full threshold scan and then choose the most appropriate range for plotting and interpretation.
Input: ABGD_conspecificity_matrix.csv
Output: ABGD_threshold_scan/
Date: 2026-03-20
"""

import pandas as pd
import networkx as nx
import os

# Load matrix
matrix_path = "ABGD_conspecificity_matrix/ABGD_conspecificity_matrix.csv"
df = pd.read_csv(matrix_path, index_col=0)

# Use all strains automatically
strains = df.index.tolist()

# Output folder
output_dir = "ABGD_threshold_scan"
os.makedirs(output_dir, exist_ok=True)

# >= 50% of core genes
max_value = int(df.values.max())
min_threshold = int(max_value * 0.5)

# 100% of core genes
#max_value = int(df.values.max())
#min_threshold = int(df.values.min())

summary = []

print("Scanning thresholds from", min_threshold, "to", max_value)

for threshold in range(min_threshold, max_value + 1):

    G = nx.Graph()
    G.add_nodes_from(strains)

    for i, s1 in enumerate(strains):
        for j, s2 in enumerate(strains):
            if i < j and df.loc[s1, s2] >= threshold:
                G.add_edge(s1, s2)

    num_groups = nx.number_connected_components(G)
    summary.append((threshold, num_groups))

# Save summary
summary_df = pd.DataFrame(summary, columns=["Threshold", "Num_Groups"])
summary_path = os.path.join(output_dir, "threshold_summary.csv")
summary_df.to_csv(summary_path, index=False)

print("Done.")
print("Threshold scan saved to:", summary_path)
