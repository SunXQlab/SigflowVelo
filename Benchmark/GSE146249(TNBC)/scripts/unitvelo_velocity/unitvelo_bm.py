import unitvelo as utv
import scvelo as scv
import scanpy as sc
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd

# 1. 配置设置
# 如果有GPU，UniTVelo会自动使用。这里显式设置日志等级
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')

# 2. 加载数据
data_path = 'GSE146249(TNBC)/data/GSE169246_converted_for_python.h5ad'
adata = scv.read(data_path)
print(f"原始数据维度: {adata.shape}")

# 3. 数据清洗：处理缺失的spliced数据
# 逻辑：检查spliced层是否存在NaN值，或总计数为0的异常细胞（视“缺失”的具体定义而定）
# 针对您描述的5%缺失情况，建议直接过滤，以免影响动力学方程求解

if 'spliced' in adata.layers:
    spliced_mat = adata.layers['spliced']
    
    # 情况A: 如果数据是稀疏矩阵 (Sparse Matrix)
    if hasattr(spliced_mat, 'data'):
        # 检查是否有NaN
        has_nan = np.isnan(spliced_mat.data).any()
        if has_nan:
            # 找到包含NaN的行(细胞)
            # 注意：稀疏矩阵直接操作较复杂，这里提供一种通用过滤思路
            # 先转为稠密检测（如果内存允许），或者检查sums
            print("检测到稀疏矩阵中存在NaN，正在过滤...")
            # 简单方法：过滤掉spliced总数为0或者NaN的细胞
            sc.pp.filter_cells(adata, min_counts=1, layer='spliced')
            
    # 情况B: 如果数据是稠密矩阵 (Dense Array)
    else:
        # 检查每一行(细胞)是否存在NaN
        cell_mask = ~np.isnan(spliced_mat).any(axis=1)
        # 同时也过滤掉全0的行（无spliced表达）
        zero_mask = spliced_mat.sum(axis=1) > 0
        
        valid_cells = cell_mask & zero_mask
        
        # 执行过滤
        adata = adata[valid_cells].copy()
        print(f"过滤后数据维度: {adata.shape} (移除了不合规细胞)")
else:
    raise ValueError("未在adata中找到'spliced'层，请检查数据格式。")

# 4. 预处理 (标准 scVelo/UniTVelo 流程)
# 筛选高变基因并标准化。UniTVelo 建议保留较多基因以获得更好的拟合
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# 5. 配置 UniTVelo 模型
config = utv.config.Configuration()

# 关键参数调整
config.R2_ADJUST = True  # 使用调整后的R2，通常效果更好
config.FIT_OPTION = '1'  # '1' = Unified-time mode (统一时间模式)
                         # 这种模式非常适合分析单一连续轨迹（如敏感->耐药的转变）
# config.GPU = 0         # 默认使用0号GPU，如需指定可取消注释

# 6. 运行模型
# 注意：'label' 参数用于初始化和绘图时的着色。
# 请将 'cell_type' 替换为您数据中实际表示细胞类型/状态的列名（例如 'condition', 'leiden', 'clusters' 等）
label_key = 'condition' 

# 检查列是否存在，不存在则使用默认聚类
if label_key not in adata.obs.columns:
    print(f"Warning: '{label_key}' not found. Running simple clustering...")
    sc.tl.leiden(adata)
    label_key = 'leiden'

print(f"开始运行 UniTVelo，使用标签列: {label_key}")
adata = utv.run_model(adata, label_key, config_file=config)

print("正在检查并修复 obsm 数据格式...")

# 1. 强制修复 UMAP
if 'UMAP' in adata.obsm:
    if isinstance(adata.obsm['UMAP'], pd.DataFrame):
        print("发现 UMAP 为 DataFrame，正在转换为 Numpy Array...")
        adata.obsm['UMAP'] = adata.obsm['UMAP'].values
    # 双重保险：有时候是类似 DataFrame 的结构
    elif not isinstance(adata.obsm['UMAP'], np.ndarray):
        adata.obsm['UMAP'] = np.array(adata.obsm['UMAP'])

# 2. 通用修复：检查所有 obsm 并转换
for key in list(adata.obsm.keys()):
    if isinstance(adata.obsm[key], pd.DataFrame):
        print(f"转换 {key} 为 Numpy Array")
        adata.obsm[key] = adata.obsm[key].values
# --- 修复代码结束 ---
save_dir = 'GSE146249(TNBC)/scripts/unitvelo_velocity/figures'
# 再运行绘图代码
plt.figure(figsize=(8, 6))
scv.pl.velocity_embedding_stream(
    adata, 
    basis='UMAP', 
    color=label_key, 
    title='UniTVelo: Sensitive to Resistant Transition',
    show=False
)
plt.savefig(os.path.join(save_dir, 'unitvelo_stream.png'), dpi=300, bbox_inches='tight')
plt.close()
# 7.2 潜伏时间 (Latent Time) 分布
# UniTVelo 计算出的统一潜伏时间，能很好地量化“敏感->耐药”的进程
plt.figure(figsize=(8, 6))
scv.pl.scatter(
    adata, 
    color='latent_time', 
    color_map='gnuplot', 
    basis='UMAP',
    size=80, 
    title='UniTVelo Latent Time (Pseudotime)',
    dpi=120,
    show=False
)

plt.savefig(os.path.join(save_dir, 'unitvelo_latent_time.png'), dpi=300, bbox_inches='tight')
plt.close() # 关闭画布释放内存

print(f"图片已保存至: {save_dir}")
# 8. 保存结果
output_path = 'GSE146249(TNBC)/scripts/unitvelo_velocity/MonoMacro_Part1_with_velocity.h5ad'
adata.write(output_path, compression='gzip')
print(f"分析完成，结果已保存至: {output_path}")