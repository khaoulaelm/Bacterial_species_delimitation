import pandas as pd
import networkx as nx
import os

# Load matrix
df = pd.read_csv("conspecificity_matrix.csv", index_col=0)

# Use all strains automatically
strains = df.index.tolist()

# Output folder
os.makedirs("threshold_scan", exist_ok=True)

# >= 50% of core genes
max_value = int(df.values.max())
min_threshold = int(max_value * 0.5)

summary = []

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
summary_df.to_csv("threshold_scan/threshold_summary.csv", index=False)

print("Done. Threshold scan saved.")
