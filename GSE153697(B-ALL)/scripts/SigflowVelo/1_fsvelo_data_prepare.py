import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv
import torch
import scipy
from scipy.spatial import distance_matrix


adata = sc.read('../../data/scPDAC_Malignant_cells_final.h5ad')

gene_node_df = adata.uns['flowsig_network']['flow_var_info']['Type']

# count矩阵
flow_vars = adata.uns['flowsig_network']['network']['flow_vars']
x_flow_df = pd.DataFrame(adata.obsm['X_flow'], index=adata.obs.index, columns=flow_vars).T
(num_outflow,num_inflow,num_gem)=(5,23,5)
x_flow_df.index = (
    [f'Outflow {i+1}' for i in range(num_outflow)] +
    [f'Inflow {i+1}' for i in range(num_inflow)] +
    list(x_flow_df.index[num_outflow+num_inflow:])
)
x_flow_df = x_flow_df.T
print(x_flow_df.columns)
x_flow_df.to_csv('../../resluts/SigflowVelo/adata/intermediate_results/count_file.csv')

# 查看网络
adjacency = adata.uns['flowsig_network']['network']['adjacency']
adjacency_df = pd.DataFrame(adjacency, index=x_flow_df.columns, columns=x_flow_df.columns)
print(adjacency_df)

# gem-outflow的连接关系
gem_outflow_df = adjacency_df.iloc[:num_outflow, -num_gem:]
gem_outflow_matrix = adjacency_df.iloc[:num_outflow, -num_gem:]
print(gem_outflow_matrix)
gem_outflow_matrix.to_csv('../../resluts/SigflowVelo/adata/Intermediate_results/Gem-Outflow_score_flie.csv')

# Gem-Outflow_Link_file
connections = []
for gem in gem_outflow_df.columns:
    outflows = gem_outflow_df.index.tolist()
    for outflow in outflows:
        connections.append({'GEM': gem, 'Outflow': outflow})
connections_df = pd.DataFrame(connections)

gem_outflow_df = connections_df
print(gem_outflow_df)
gem_outflow_df.to_csv('../../resluts/SigflowVelo/adata/Intermediate_results/Gem-Outflow_link_flie.csv')

# Inflow-Gem_Link_file

inflow_gem_df = adjacency_df.iloc[-num_gem:,num_outflow:num_outflow+num_inflow]

connections = []
for gem in inflow_gem_df.columns:
    outflows = inflow_gem_df.index[inflow_gem_df[gem] > 0].tolist()
    for outflow in outflows:
        connections.append({'Inflow': gem, 'Gem': outflow})
connections_df = pd.DataFrame(connections)

inflow_gem_df = connections_df
print(inflow_gem_df)
inflow_gem_df.to_csv('../../resluts/SigflowVelo/adata/Intermediate_results/Inflow-Gem_link_flie.csv')

# inflow-gem的连接信号强度
inflow_gem_df = adjacency_df.iloc[-num_gem:,num_outflow:num_outflow+num_inflow]
print(inflow_gem_df)
inflow_gem_df.to_csv('../../resluts/SigflowVelo/adata/Intermediate_results/Inflow-Gem_score_flie.csv')

