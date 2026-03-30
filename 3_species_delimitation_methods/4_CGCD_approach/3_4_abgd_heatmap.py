#!/usr/bin/env python3
"""
Author: Khaoula El Mchachti
Description: Plot a clustered heatmap (clustermap) of the ABGD conspecificity matrix.
Input: ABGD_conspecificity_matrix.csv
Output: ABGD_heatmap.pdf
Date: 2026-03-20
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# Input and output paths
matrix_path = os.path.expanduser("~/Bacterial_species_delimitation/3_species_delimitation_methods/4_CGCD_approach/ABGD_conspecificity_matrix/ABGD_conspecificity_matrix.csv")

# Directory where the file will be saved
output_dir = os.path.expanduser("~/Bacterial_species_delimitation/3_species_delimitation_methods/4_CGCD_approach/ABGD_plots")

# Create the directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# File path inside the directory
output_path = os.path.join(output_dir, "ABGD_heatmap.pdf")

# Load matrix
df = pd.read_csv(matrix_path, index_col=0)

# Set up figure size
#plt.figure(figsize=(20, 20))

# Generate heatmap with clustering
print("Generating clustered heatmap...")
g = sns.clustermap(df, 
               cmap="coolwarm",  # Better contrast
               linewidths=0.2,   # Add gridlines
               linecolor="black",# Grid color
               cbar_kws={"shrink": 0.5},  # Smaller colorbar
               xticklabels=True,  # Show x labels
               yticklabels=True,  # Show y labels
               figsize=(20, 20),  # Size
               annot=False,       # Set to True if you want numbers in cells
               fmt=".0f",         # No decimals for annotations
               dendrogram_ratio=(0.1, 0.1) # Reduce dendrogram size
              )
        

cbar_pos = g.cax.get_position()
g.cax.set_position([
    cbar_pos.x0,          
    cbar_pos.y0 + 0.1,         
    cbar_pos.width * 0.5,  
    cbar_pos.height *0.5     
])

print("Saving heatmap to:", output_path)
g.savefig(output_path, format='pdf')
plt.close()
