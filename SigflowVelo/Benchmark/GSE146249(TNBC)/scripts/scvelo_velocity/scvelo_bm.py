import scvelo as scv
import scanpy as sc
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd

# 1. 配置
# 自动检测CPU核心数并设定合理的并行数 (例如使用 75% 的核心)
total_cores = os.cpu_count()
n_jobs = max(4, int(total_cores * 0.5)) 
print(f"检测到 {total_cores} 个核心，本次 scVelo 分析将使用 {n_jobs} 个核心。")
# 调高 verbosity 以便看到拟合进度
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')

# 设置保存路径
save_dir = 'GSE146249(TNBC)/scripts/scvelo_velocity/figures'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# 2. 加载数据
data_path = 'GSE146249(TNBC)/data/GSE169246_converted_for_python.h5ad'
adata = scv.read(data_path)
print(f"原始维度: {adata.shape}")

# 3. 数据清洗 (必须与 UnitVelo 保持一致)
print(f"清洗前维度: {adata.shape}")
if 'spliced' in adata.layers:
    spliced_mat = adata.layers['spliced']
    
    # 1. 计算每个细胞的 spliced 总数 (兼容稀疏和稠密矩阵)
    if hasattr(spliced_mat, 'sum'):
        spliced_counts = spliced_mat.sum(axis=1)
        # 如果是 matrix 对象 (稀疏矩阵常见)，转为 flat array
        if hasattr(spliced_counts, 'A1'):
            spliced_counts = spliced_counts.A1
        else:
            spliced_counts = np.array(spliced_counts).flatten()
    else:
        # 稠密矩阵
        spliced_counts = np.sum(spliced_mat, axis=1)

    # 2. 生成过滤掩码 (Spliced 计数 > 0)
    # 也可以处理可能的 NaN
    spliced_counts = np.nan_to_num(spliced_counts)
    cell_mask = spliced_counts > 0
    
    # 3. 执行过滤
    adata = adata[cell_mask].copy()
    print(f"清洗后维度: {adata.shape} (已移除 {np.sum(~cell_mask)} 个无 spliced 数据的细胞)")
    
else:
    print("Warning: 未找到 'spliced' 层，跳过该层过滤。")
# 4. 预处理
# 筛选高变基因 (scVelo 动力学模式通常需要较严格的筛选)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# 5. 运行动力学模型 (Dynamical Mode)
# 这一步比 UnitVelo 慢，因为要逐个基因解微分方程
print("正在恢复动力学参数 (Recovering Dynamics)...")
scv.tl.recover_dynamics(adata, n_jobs=8) # 根据服务器核数调整 n_jobs

print("计算 Velocity...")
scv.tl.velocity(adata, mode='dynamical')

print("计算 Velocity Graph...")
scv.tl.velocity_graph(adata)

# 6. 计算 Latent Time (潜伏时间)
# 这是 scVelo 动力学模式的核心优势，能精准量化细胞在轨迹上的位置
print("计算 Latent Time...")
scv.tl.latent_time(adata)
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
# 7. 绘图与保存
label_key = 'condition' # 请替换为你的真实列名，如 'cell_type', 'leiden'
if label_key not in adata.obs.columns:
    sc.tl.leiden(adata)
    label_key = 'leiden'

# 7.1 流线图
scv.pl.velocity_embedding_stream(
    adata, basis='UMAP', color=label_key, 
    title='scVelo (Dynamical): sensitive -> resistant',
    show=False
)
plt.savefig(os.path.join(save_dir, 'scvelo_stream.png'), dpi=300, bbox_inches='tight')
plt.close()

# 7.2 Latent Time
scv.pl.scatter(
    adata, color='latent_time', color_map='gnuplot', size=80,
    title='scVelo Latent Time',basis='UMAP',
    show=False
)
plt.savefig(os.path.join(save_dir, 'scvelo_latent_time.png'), dpi=300, bbox_inches='tight')
plt.close()

# 8. 保存结果
out_path = 'GSE146249(TNBC)/scripts/scvelo_velocity/MonoMacro_Part1_scvelo.h5ad'
adata.write(out_path, compression='gzip')
print(f"scVelo 分析完成，已保存至 {out_path}")