import pandas as pd
import numpy as np
import os
import scanpy as sc
import scvelo as scv
import torch
import scipy
from scipy.spatial import distance_matrix


#加载数据
print('Loading files...')
adata = sc.read('result\script_2_result\Mono-Macro_PacTissue_immune.h5ad')
print('Loading files complete.')

#变量准备
## 导入adata.uns['flowsig_network']['flow_var_info']信息
try:
    gene_node_df = adata.uns['flowsig_network']['flow_var_info']['Type']
except KeyError:
    print("错误：在 anndata 文件的 adata.uns['flowsig_network']['flow_var_info'] 中未找到 'Type'。")
    print("请确保您加载的 .h5ad 文件是 脚本 2 (2_flowsig.py) 的最终、完整输出。")
    raise

## 使用 .value_counts() 自动统计每种类型的确切数量
variable_counts = gene_node_df.value_counts()
# 从统计结果中安全地获取每种类型的数量
# 使用 .get(key, 0) 确保即使某个类型不存在（数量为0），脚本也不会报错
outflow_count = variable_counts.get('outflow', 0)
inflow_count = variable_counts.get('inflow', 0)
module_count = variable_counts.get('module', 0) # 假设 脚本 2 将模块命名为 'GEM'

# 打印出来进行验证
print("--- 动态检测到网络变量 ---")
print(f"Outflow 数量 (outflow_count): {outflow_count}")
print(f"Inflow 数量 (inflow_count): {inflow_count}")
print(f"GEM 数量 (module_count): {module_count}")
total_vars = outflow_count + inflow_count + module_count
print(f"检测到的总变量数: {total_vars}")
print(f"adata中的总变量数: {len(gene_node_df)}")
print("------------------------------")

# 检查总数是否匹配
if total_vars != len(gene_node_df):
    print(f"警告：统计到的变量总数 ({total_vars}) 与 AnnData 中的变量总数 ({len(gene_node_df)}) 不匹配！")
    print("这可能是因为 'Type' 列中存在未知的变量类型。")
    print("观测到的所有类型:", variable_counts.index.tolist())

# 提取count矩阵，从adata.obsm['X_flow'] 提取了所有网络变量（Outflows, Inflows, GEMs）在每个细胞中的“表达”或“活性”水平。
flow_vars = adata.uns['flowsig_network']['network']['flow_vars']
x_flow_df = pd.DataFrame(adata.obsm['X_flow'], index=adata.obs.index, columns=flow_vars).T



# 修改行名：这里要改成你自己的数据集的Outflow和Inflow和Gem的个数
x_flow_df.index = (
    [f'Outflow {i+1}' for i in range(outflow_count)] + 
    [f'Inflow {i+1}' for i in range(inflow_count)] + 
    list(x_flow_df.index[outflow_count + inflow_count:])
)
x_flow_df = x_flow_df.T
print(x_flow_df.columns)
x_flow_df.to_csv('result/script_3_result/Mono-Macro_count_file.csv')

# 提取邻接矩阵，这是一个方阵，行和列都是所有的网络变量(Outflows, Inflows, GEMs),值代表链接强度。
print('Extracting adjacency matrix...')
adjacency = adata.uns['flowsig_network']['network']['adjacency']
adjacency_df = pd.DataFrame(adjacency, index=x_flow_df.columns, columns=x_flow_df.columns)
print('Extraction complete.')
print(adjacency_df)

# 拆分网络
## 切片邻接矩阵，代表了从GEMs到Outflows的连接。
print('Sclicing connection from GEMs to Outflows...')
gem_outflow_df = adjacency_df.iloc[:outflow_count, -module_count:]
gem_outflow_matrix = adjacency_df.iloc[:outflow_count, -module_count:]
print(gem_outflow_matrix)
gem_outflow_matrix.to_csv('result/script_3_result/Mono-Macro_Gem-Outflow_score_flie.csv')

## Gem-Outflow_Link_file
connections = []
for gem in gem_outflow_df.columns:
    outflows = gem_outflow_df.index[gem_outflow_df[gem] > 0].tolist()
    for outflow in outflows:
        connections.append({'GEM': gem, 'Outflow': outflow})
connections_df = pd.DataFrame(connections)

gem_outflow_df = connections_df
print(gem_outflow_df)
gem_outflow_df.to_csv('result/script_3_result/Mono-Macro_Gem-Outflow_link_flie.csv')

# 切片邻接矩阵，代表了从Inflows到GEMs的连接。
inflow_gem_df = adjacency_df.iloc[-module_count:,outflow_count:outflow_count+inflow_count] #行是GEMs，列是inflow
## Inflow-Gem_Link_file
connections = []
for inflow in inflow_gem_df.columns:
    gems = inflow_gem_df.index[inflow_gem_df[inflow] > 0].tolist()
    for gem in gems:
        connections.append({'Inflow': inflow, 'Gem': gem})
connections_df = pd.DataFrame(connections)

inflow_gem_df_link = connections_df
print(inflow_gem_df_link)
inflow_gem_df_link.to_csv('result/script_3_result/Mono-Macro_Inflow-Gem_link_flie.csv')

# inflow-gem的连接信号强度
inflow_gem_df = adjacency_df.iloc[-module_count:,outflow_count:outflow_count+inflow_count]
print(inflow_gem_df)
inflow_gem_df.to_csv('result/script_3_result/Mono-Macro_Inflow-Gem_score_flie.csv')

