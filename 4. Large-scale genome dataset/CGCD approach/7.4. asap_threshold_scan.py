import pandas as pd
import networkx as nx
import os

# === Load your files ===
matrix_path = "~/asap_conspecificity_matrix.csv"
vub_group_path = "~/vub_groups.csv"

df = pd.read_csv(os.path.expanduser(matrix_path), index_col=0)
vub_df = pd.read_csv(os.path.expanduser(vub_group_path))
vub_strains = vub_df["Strain"].tolist()

# Create submatrix for VUB strains
df_vub = df.loc[vub_strains, vub_strains]

# === Create output directory ===
output_dir = "ASAP_groupings"
os.makedirs(output_dir, exist_ok=True)

# === Loop through thresholds ===
group_summary = []

for threshold in range(178, 223): #(>= 80%)
    G = nx.Graph()
    G.add_nodes_from(df.index)

    # Add edges if pair is supported by ≥ threshold genes
    for i, strain_i in enumerate(df.index):
        for j, strain_j in enumerate(df.columns):
            if i < j and df.loc[strain_i, strain_j] >= threshold:
                G.add_edge(strain_i, strain_j)

    # Extract connected components
    components = list(nx.connected_components(G))
    
    total_num_groups = len(components)


    # Track number of VUB groups only
    vub_groups = [comp for comp in components if any(s in vub_strains for s in comp)]
    vub_only_components = []

    for comp in vub_groups:
        vub_subgroup = set(comp).intersection(set(vub_strains))
        if vub_subgroup:
            vub_only_components.append(vub_subgroup)

    num_vub_groups = len(vub_only_components)
    group_summary.append((threshold, num_vub_groups, total_num_groups))


    # === Save if exactly 11 VUB groups found ===
    if num_vub_groups == 11:
        # Save VUB-only groupings
        vub_group_assignments = []
        for group_id, members in enumerate(vub_only_components, start=1):
            for strain in members:
                vub_group_assignments.append({"Strain": strain, "Group": group_id})

        vub_group_df = pd.DataFrame(vub_group_assignments)
        vub_out_path = os.path.join(output_dir, f"vub_ASAPgroups_t{threshold}.csv")
        vub_group_df.to_csv(vub_out_path, index=False)

        # Save full groupings
        full_group_assignments = []
        for group_id, members in enumerate(components, start=1):
            for strain in members:
                full_group_assignments.append({"Strain": strain, "Group": group_id})

        full_group_df = pd.DataFrame(full_group_assignments)
        full_out_path = os.path.join(output_dir, f"all_strain_ASAPgroups_t{threshold}.csv")
        full_group_df.to_csv(full_out_path, index=False)

# === Save summary ===
summary_df = pd.DataFrame(group_summary, columns=["Threshold", "Num_VUB_Groups", "Num_Total_Groups"])
summary_df.to_csv(os.path.join(output_dir, "ASAP_groupings_summary.csv"), index=False)

print(f"All groupings and summary saved to: {output_dir}")
