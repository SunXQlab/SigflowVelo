import pandas as pd
import numpy as np
import os
import scanpy as sc
import scvelo as scv
import torch
import scipy
from scipy.spatial import distance_matrix

# Extract network information from FlowSig analysis results
# Prepare structured data for subsequent flow velocity inference model

# Read the final analysis results from the previous step
adata = sc.read("SigflowVelo/benchmark/GSE164551(MM)/data/SigflowVelo_data/NK_final.h5ad")

# Get actual flow variable count and types from FlowSig analysis results
flow_var_info = adata.uns['flowsig_network']['flow_var_info']

# Count the number of each type
outflow_count = sum(flow_var_info['Type'] == 'outflow')
inflow_count = sum(flow_var_info['Type'] == 'inflow') 
gem_count = sum(flow_var_info['Type'] == 'module')

print(f"Outflow count: {outflow_count}")
print(f"Inflow count: {inflow_count}")
print(f"GEM count: {gem_count}")

gene_node_df = adata.uns['flowsig_network']['flow_var_info']['Type']

# Count matrix
flow_vars = adata.uns['flowsig_network']['network']['flow_vars']
x_flow_df = pd.DataFrame(adata.obsm['X_flow'], index=adata.obs.index, columns=flow_vars).T

# Modify row names: dynamically generate based on actual counts
x_flow_df.index = (
    [f'Outflow {i+1}' for i in range(outflow_count)] + 
    [f'Inflow {i+1}' for i in range(inflow_count)] + 
    list(x_flow_df.index[outflow_count + inflow_count:])
)
x_flow_df = x_flow_df.T
print(x_flow_df.columns)

# Save expression/activity levels of each cell on various flow variables
# This serves as input features for deep learning models
x_flow_df.to_csv("SigflowVelo/benchmark/GSE164551(MM)/data/SigflowVelo_data/count_file.csv")

# View network
adjacency = adata.uns['flowsignetwork']['network']['adjacency']
adjacency_df = pd.DataFrame(adjacency, index=x_flow_df.columns, columns=x_flow_df.columns)
print(adjacency_df)

# GEM-outflow connections (dynamic indexing based on actual counts)
# Note: Assuming GEMs come after outflows and inflows in the adjacency matrix
gem_start_idx = outflow_count + inflow_count
gem_end_idx = gem_start_idx + gem_count
gem_outflow_df = adjacency_df.iloc[gem_start_idx:gem_end_idx, :outflow_count]
gem_outflow_matrix = adjacency_df.iloc[gem_start_idx:gem_end_idx, :outflow_count]
print(gem_outflow_matrix)

# Save GEM to Outflow connection weight matrix
gem_outflow_matrix.to_csv('SigflowVelo/benchmark/GSE164551(MM)/data/SigflowVelo_data/Gem-Outflow_score_file.csv')

# Gem-Outflow_Link_file
connections = []
for gem in gem_outflow_df.columns:
    outflows = gem_outflow_df.index[gem_outflow_df[gem] > 0].tolist()
    for outflow in outflows:
        connections.append({'GEM': gem, 'Outflow': outflow})
connections_df = pd.DataFrame(connections)

gem_outflow_link_df = connections_df
print(gem_outflow_link_df)

# Save effective GEM to Outflow connection list
gem_outflow_link_df.to_csv('SigflowVelo/benchmark/GSE164551(MM)/data/SigflowVelo_data/Gem-Outflow_link_file.csv')

# Inflow-Gem_Link_file
inflow_gem_df = adjacency_df.iloc[:inflow_count, gem_start_idx:gem_end_idx]

connections = []
for gem in inflow_gem_df.columns:
    inflows = inflow_gem_df.index[inflow_gem_df[gem] > 0].tolist()
    for inflow in inflows:
        connections.append({'Inflow': inflow, 'Gem': gem})
connections_df = pd.DataFrame(connections)

inflow_gem_link_df = connections_df
print(inflow_gem_link_df)

# Save effective Inflow to GEM connection list
inflow_gem_link_df.to_csv('SigflowVelo/benchmark/GSE164551(MM)/data/SigflowVelo_data/Inflow-Gem_link_file.csv')

# Inflow-Gem connection signal strength
inflow_gem_score_df = adjacency_df.iloc[:inflow_count, gem_start_idx:gem_end_idx]
print(inflow_gem_score_df)

# Save Inflow to GEM connection weight matrix
inflow_gem_score_df.to_csv('SigflowVelo/benchmark/GSE164551(MM)/data/SigflowVelo_data/Inflow-Gem_score_file.csv')