import pandas as pd
import networkx as nx
import os

# Load conspecificity matrix
df = pd.read_csv("conspecificity_matrix.csv", index_col=0)

# All strains
strains = df.index.tolist()

# Plateau range (from your plot)
START = 1494
END = 2026

# Output directory
os.makedirs("groups_plateau", exist_ok=True)

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
        f"groups_plateau/groups_t{threshold}.csv",
        index=False
    )

print("Done. Groups extracted for all thresholds in the plateau.")
