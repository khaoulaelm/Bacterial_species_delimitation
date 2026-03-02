#!/usr/bin/env python3

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
from matplotlib import rcParams

# Load your data
df = pd.read_csv("fastani_output.csv", sep="\t", header=None)
df.columns = ["genome1", "genome2", "ani", "fragments", "total"]

# Extract just the strain name (remove full path and extension)
df["genome1"] = df["genome1"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
df["genome2"] = df["genome2"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])

# Pivot to wide matrix format
ani_matrix = df.pivot(index="genome1", columns="genome2", values="ani")

# Fill missing values with 0 
ani_matrix = ani_matrix.fillna(0)

# Ensure all strain names are present in both rows and columns
all_strains = sorted(set(ani_matrix.index).union(set(ani_matrix.columns)))
ani_matrix = ani_matrix.reindex(index=all_strains, columns=all_strains)

# Force symmetry
ani_matrix = (ani_matrix + ani_matrix.T) / 2

# Convert ANI to distance
distance_matrix = 100 - ani_matrix

# Create condensed distance matrix for clustering
condensed_dist = squareform(distance_matrix)
linkage_matrix = linkage(condensed_dist, method="average")

# Plot heatmap with clustering
sns.set(style="white")
rcParams['figure.figsize'] = 25, 25

clustermap = sns.clustermap(
    ani_matrix,
    row_linkage=linkage_matrix,
    col_linkage=linkage_matrix,
    cmap="Spectral_r",        
    linewidths=0.2,
    linecolor="black",
    figsize=(25, 25),
    xticklabels=True,
    yticklabels=True,
    cbar_kws={"label": "ANI (%)"},
)

# Rotate labels for clarity
plt.setp(clustermap.ax_heatmap.get_xticklabels(), rotation=90)
plt.setp(clustermap.ax_heatmap.get_yticklabels(), rotation=0)

# Save the heatmap
output_file = "ANI_heatmap.pdf"
plt.savefig(output_file, format='pdf', bbox_inches='tight')
plt.close()
