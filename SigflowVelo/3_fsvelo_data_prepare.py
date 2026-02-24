import pandas as pd
import numpy as np
import os
import scanpy as sc
import scvelo as scv
import torch
import scipy
from scipy.spatial import distance_matrix


adata = sc.read('../output/dutcal/scPDAC_dutcal_final.h5ad')

gene_node_df = adata.uns['flowsig_network']['flow_var_info']['Type']

# count矩阵
flow_vars = adata.uns['flowsig_network']['network']['flow_vars']
x_flow_df = pd.DataFrame(adata.obsm['X_flow'], index=adata.obs.index, columns=flow_vars).T
# 修改行名：这里要改成你自己的数据集的Outflow和Inflow和Gem的个数
x_flow_df.index = (
    [f'Outflow {i+1}' for i in range(27)] + 
    [f'Inflow {i+1}' for i in range(31)] + 
    list(x_flow_df.index[58:])
)
x_flow_df = x_flow_df.T
print(x_flow_df.columns)
x_flow_df.to_csv('../CCCInputData/Dutcal/count_file.csv')

# 查看网络
adjacency = adata.uns['flowsig_network']['network']['adjacency']
adjacency_df = pd.DataFrame(adjacency, index=x_flow_df.columns, columns=x_flow_df.columns)
print(adjacency_df)

# gem-outflow的连接关系
gem_outflow_df = adjacency_df.iloc[:20, -10:]
gem_outflow_matrix = adjacency_df.iloc[:20, -10:]
print(gem_outflow_matrix)
gem_outflow_matrix.to_csv('../CCCInputData/Dutcal/Gem-Outflow_score_flie.csv')

# Gem-Outflow_Link_file
connections = []
for gem in gem_outflow_df.columns:
    outflows = gem_outflow_df.index[gem_outflow_df[gem] > 0].tolist()
    for outflow in outflows:
        connections.append({'GEM': gem, 'Outflow': outflow})
connections_df = pd.DataFrame(connections)

gem_outflow_df = connections_df
print(gem_outflow_df)
gem_outflow_df.to_csv('../CCCInputData/Dutcal/Gem-Outflow_link_flie.csv')

# Inflow-Gem_Link_file

inflow_gem_df = adjacency_df.iloc[-10:,27:58]

connections = []
for gem in inflow_gem_df.columns:
    outflows = inflow_gem_df.index[inflow_gem_df[gem] > 0].tolist()
    for outflow in outflows:
        connections.append({'Inflow': gem, 'Gem': outflow})
connections_df = pd.DataFrame(connections)

inflow_gem_df = connections_df
print(inflow_gem_df)
inflow_gem_df.to_csv('../CCCInputData/Dutcal/Inflow-Gem_link_flie.csv')

# inflow-gem的连接信号强度
inflow_gem_df = adjacency_df.iloc[-10:,27:58]
print(inflow_gem_df)
inflow_gem_df.to_csv('../CCCInputData/Dutcal/Inflow-Gem_score_flie.csv')

