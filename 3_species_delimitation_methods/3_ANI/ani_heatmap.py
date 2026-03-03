#!/usr/bin/env python3

"""
Author: Khaoula El Mchachti
Description: Generate clustered ANI heatmap from FastANI output
Input: fastani_results.csv 
Output: ANI_heatmap.pdf
Date: 2026-03-02
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
from matplotlib import rcParams

# Load FastANI output
input_file = "ani_results/fastani_results.csv"
output_file = "ani_results/ANI_heatmap.pdf"

print("Loading FastANI results...")


df = pd.read_csv(input_file, sep="\t", header=None)
df.columns = ["genome1", "genome2", "ani", "fragments", "total"]

# Extract just the strain name (remove full path and extension)
df["genome1"] = df["genome1"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
df["genome2"] = df["genome2"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])

# Create ANI matrix
ani_matrix = df.pivot(index="genome1", columns="genome2", values="ani")
ani_matrix = ani_matrix.fillna(0)

# Ensure symmetry and full matrix
all_strains = sorted(set(ani_matrix.index).union(set(ani_matrix.columns)))
ani_matrix = ani_matrix.reindex(index=all_strains, columns=all_strains)
ani_matrix = (ani_matrix + ani_matrix.T) / 2

# Clustering
print("Performing hierarchical clustering...")

distance_matrix = 100 - ani_matrix
condensed_dist = squareform(distance_matrix) # Create condensed distance matrix for clustering
linkage_matrix = linkage(condensed_dist, method="average")

# Plot heatmap 

print("Generating heatmap...")

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

plt.savefig(output_file, format='pdf', bbox_inches='tight')
plt.close()

print(f"Heatmap saved to: {output_file}")
