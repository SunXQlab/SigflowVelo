import scanpy as sc
import pandas as pd
import flowsig as fs
import numpy as np

#变量准备
condition_key = 'condition'
#创建新版本函数需要的实例
my_flow_config = fs.pp.FlowSigConfig(
    gem_expr_key='X_gem',           
    scale_gem_expr=False,          
    flowsig_network_key='flowsig_network', 
    flowsig_expr_key='X_flow'       
)

#加载文件
##.rds直接转换的.h5ad文件，包含全部信息。
print('Loading .h5ad file')
adata = sc.read('GSE146249(TNBC)\data\GSE169246_converted_for_python.h5ad')
# adata = sc.read('result/PacTissue_immune+GEMs.h5ad')
print('File loading complete.')

##由上一个脚本给出的细胞通讯分析。
print('Loading .csv file')
cellchat_resistant= pd.read_csv('GSE146249(TNBC)\data\PacTissue_celltype_communications_resistant.csv')
cellchat_sensitive= pd.read_csv('GSE146249(TNBC)\data\PacTissue_celltype_communications_sensitive.csv')
print('File loading complete.')


#筛选交互。
print("---正在筛选 Mono-Macro 细胞之间的交互 ---")

## 首先检查'celltype'列是否存在
if 'celltype' not in adata.obs.columns:
    print(f"错误: 在 .h5ad 文件的元数据 (adata.obs) 中未找到 'celltype' 列。")
    print(f"R 脚本 1_S2Rsig.R 指定了 group.by = 'celltype'。请确保此列已正确导出。")
    print(f"adata.obs 中可用的列: {adata.obs.columns.tolist()}")
    raise KeyError("未在 adata.obs 中找到 'celltype' 列。")

## 获取 .h5ad 文件中所有的细胞类型名称
all_celltypes = adata.obs['celltype'].unique().tolist()
print(f"在 .h5ad 文件中找到的所有细胞类型: {all_celltypes}")

## 定义T细胞组和巨噬细胞组。
t_celltype = ['Mono-Macro']
macro_celltype = ['Mono-Macro']

## 筛选 "T -> M" 或 "M -> T" 的交互
cellchat_sensitive_filtered = cellchat_sensitive[
    (cellchat_sensitive['source'].isin(t_celltype) & cellchat_sensitive['target'].isin(macro_celltype)) |
    (cellchat_sensitive['source'].isin(macro_celltype) & cellchat_sensitive['target'].isin(t_celltype))
]
cellchat_resistant_filtered = cellchat_resistant[
    (cellchat_resistant['source'].isin(t_celltype) & cellchat_resistant['target'].isin(macro_celltype)) |
    (cellchat_resistant['source'].isin(macro_celltype) & cellchat_resistant['target'].isin(t_celltype))
]
print(f"筛选 'sensitive' 后保留的交互数: {len(cellchat_sensitive)}")
print(f"筛选 'resistant' 后保留的交互数: {len(cellchat_resistant)}")
### 将筛选后的结果保存到新文件，方便检查
print('Saving filtered interactions...')
cellchat_sensitive_filtered.to_csv('result/script_2_result/Mono-Macro_cellchat_sensitive.csv', index=False)
cellchat_resistant_filtered.to_csv('result/script_2_result/Mono-Macro_cellchat_resistant.csv', index=False)
## 确保这些的键与它们的条件标签对齐
cellchat_output_key = 'cellchat_output'


## 不对作用进行筛选，直接纳入计算。
adata.uns[cellchat_output_key] = {'sensitive': cellchat_sensitive,
                                  'resistant': cellchat_resistant}
print('Saving interactions complete.')



#把基因表达矩阵添加到layers中
print('Prepare "counts" layers and index...')
counts_matrix = adata.X.copy() 
adata.layers["counts"] = counts_matrix

# 设置细胞索引的名称
if adata.obs.index.name is None:
    adata.obs.index.name = 'cell_id' 

# 设置基因索引的名称
if adata.var.index.name is None:
    adata.var.index.name = 'gene_id' 

# 使用 pyLIGER 从未归一化的基因表达计数中构建了 10 个基因表达模块（GEMs）
print('Construct GEMs by pyLIGER')
fs.pp.construct_gems_using_pyliger(adata,
                                n_gems = 10,
                                layer_key = 'counts',
                                condition_key = condition_key)

gems_result = adata.obsm['X_gem']
gems_df = pd.DataFrame(gems_result, index=adata.obs_names, columns=[f"Module_{i+1}" for i in range(gems_result.shape[1])])
gems_df.to_csv("result/script_2_result/Mono-Macro_gems_result.csv")
adata.write_h5ad("result/script_2_result/PacTissue_immune_Mono-Macro_GEMs.h5ad")  


print('Setup FlowSig...')
#核心步骤一，目的是将“基因表达”转换为“信号流表达”:配体-受体-转录因子
print('Construct signal flow from CellChat...')
fs.pp.construct_flows_from_cellchat(adata,
                                cellchat_output_key,
                                model_organism = 'human',
                                config = my_flow_config ,
                                )

flow_key = my_flow_config.flowsig_expr_key  #'X_flow'：n_cells * n_flows的一个矩阵
flow_data = adata.obsm[flow_key]
## 执行 log(x+1)
adata.obsm[flow_key] = np.log1p(flow_data)

# 测试：检查总共有多少星系变量
# --- 新增：检查总变量池大小 ---
n_gems = adata.obsm['X_gem'].shape[1]
n_flows = adata.obsm['X_flow'].shape[1]
total_vars = n_gems + n_flows

print("--- 检查：总变量池（筛选前）---")
print(f"基因模块 (GEMs) 数量: {n_gems}")
print(f"信号流 (Flows) 数量: {n_flows}")
print(f"总变量 (可供筛选的) 数量: {total_vars}")
print("---------------------------------")


#核心步骤二，是一个特征筛选步骤，只保留敏感/抗药之间存在显著差异的变量（无论是GEMs还是flows）。
print('Filtering GEMs and Flows by apparent difference between sensitive group and resistant group...')
fs.pp.determine_informative_variables(adata,  
                                    config = my_flow_config,
                                    spatial = False,
                                    condition_key = condition_key,
                                    control = 'sensitive',
                                    qval_threshold = 0.05,
                                    logfc_threshold = 0.5,
                                    )


# 测试：检查一下我们找到了多少变量
flow_var_info = adata.uns['flowsig_network']['flow_var_info']
print(f"--- 成功找到了 {len(flow_var_info)} 个信息性变量 ---")

# 测试：检查一下sensitive组中有多少个细胞
# --- 新增：验证 'sensitive' 组的细胞总数 (N) ---
print("--- 正在验证 Control 组的细胞总数 (N) ---")

# condition_key 在脚本开头定义为 'condition'
# control 在 learn_intercellular_flows 中定义为 'sensitive'
n_sensitive_cells = (adata.obs[condition_key] == 'sensitive').sum()
n_resistant_cells = (adata.obs[condition_key] == 'resistant').sum()

# n_variables_p 就是你上一行打印的 变量总数
n_variables_p = len(flow_var_info) 

print(f"Control 组 ('sensitive') 细胞数 (N): {n_sensitive_cells}")
print(f"Perturbed 组 ('resistant') 细胞数: {n_resistant_cells}")
print(f"当前筛选出的信息性变量数 (P): {n_variables_p}")
print("-------------------------------------------------")

## 保存结果
flow_var_info = adata.uns['flowsig_network']['flow_var_info']
flow_var_info_df = pd.DataFrame(flow_var_info)
flow_var_info_df.to_csv('result/Mono-Macro_flow_var_info.csv', index=True)
adata.write_h5ad("result/Mono-Macro_PacTissue_immune+GEMs+Signal.h5ad")  
  

print("--- 正在清理 'X_flow' 矩阵中的 NaN 值，替换为 0.0 ---")

# np.nan_to_num 会将 NaN 替换为 0.0 (也会处理 inf)
adata.obsm[my_flow_config.flowsig_expr_key] = np.nan_to_num(
    adata.obsm[my_flow_config.flowsig_expr_key], 
    nan=0.0, 
    posinf=0.0, 
    neginf=0.0
)
print("--- NaN 清理完成 --- \n")





#核心步骤三，学习因果网络
print('Learning causation network of intercellular signal flows...')
## 使用上一步筛选出来的数据，通过自举法(Bootstrapping)和可能的结构方程模型 (SEM) 或因果发现算法来推断变量之间的有向边
fs.tl.learn_intercellular_flows(adata,
                        condition_key = condition_key,
                        control = 'sensitive', 
                        config = my_flow_config,
                        use_spatial = False,
                        n_jobs = 5,
                        n_bootstraps = 100)



##网络修正：利用先验的生物学知识（比如信号必须从L -> R -> TF）来强制修正边的方向
fs.tl.apply_biological_flow(adata,
                        flowsig_network_key = 'flowsig_network',
                        adjacency_key = 'adjacency')


##网络剪枝：移除那些“不稳定”的链接：只保留在70%以上的自举运行中都出现过的连接。
edge_threshold = 0.7
fs.tl.filter_low_confidence_edges(adata,
                                edge_threshold = edge_threshold,
                                flowsig_network_key = 'flowsig_network',
                                adjacency_key = 'adjacency')
adata.write('result/script_2_result/Mono-Macro_burkhardt21_merged.h5ad', compression='gzip')
print('Network constrcution complete.')

# 将网络保存为csv
print('Saving files...')
adjacency = adata.uns['flowsig_network']['network']['adjacency']
adjacency_df = pd.DataFrame(adjacency, index=flow_var_info_df.index, columns=flow_var_info_df.index)
adjacency_df.to_csv('result/script_2_result/Mono-Macro_adjacency_network.csv', index=True)

adata.write_h5ad("result/script_2_result/Mono-Macro_PacTissue_immune.h5ad")  