#!/usr/bin/env python3
"""
Author: Khaoula El Mchachti
Description: Extract strain groups for each threshold within a detected plateau of the conspecificity threshold scan.
Input: ASAP_conspecificity_matrix.csv
Output: ASAP_groups_plateau/groups_t{threshold}.csv
Date: 2026-03-30
Last modified: 2026-08-17
"""

import pandas as pd
import networkx as nx
import os

# Find the directory containing this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Project root directory
PROJECT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "../.."))

#Input: ASAP conspecificity matrix
matrix_file = os.path.join(
    PROJECT_DIR,
    "3_species_delimitation_methods",
    "4_CGCD_approach",
    "ASAP_conspecificity_matrix",
    "ASAP_conspecificity_matrix.csv"
)

# Load conspecificity matrix
df = pd.read_csv(matrix_file, index_col=0)

# All strains
strains = df.index.tolist()

# Plateau range obtained from the pervious threshold scan.
# Read the detected plateau
plateau_file = os.path.join(
    PROJECT_DIR,
    "3_species_delimitation_methods",
    "4_CGCD_approach",
    "ASAP_plots",
    "best_plateau.txt"
)

with open(plateau_file) as f:
   START, END = map( int, f.read().strip().split())
   
# Output directory
output_dir = os.path.join(
    PROJECT_DIR,
    "3_species_delimitation_methods",
    "4_CGCD_approach",
    "ASAP_groups_plateau"
)

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

    output_file = os.path.join(
        output_dir,
        f"groups_t{threshold}.csv"
    )

    pd.DataFrame(rows).to_csv(
        output_file,
        index=False
    )
    
print("Done. Groups extracted for all thresholds in the plateau.")
print(f"Output directory: {output_dir}")
