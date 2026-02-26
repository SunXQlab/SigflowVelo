import os
# 1. 强制设置环境变量，屏蔽 GPU (最稳妥的 CPU 模式方法)
os.environ["CUDA_VISIBLE_DEVICES"] = ""

import torch
import deepvelo as dv
import scvelo as scv
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# 设置 CPU 线程数
n_threads = max(4, int(os.cpu_count() * 0.75))
torch.set_num_threads(n_threads)
print(f"--- 正在使用 CPU 模式，线程数限制为: {n_threads} ---")

# 2. 基础路径设置
save_dir = 'GSE146249(TNBC)/scripts/deepvelo_velocity/figures'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# 3. 加载数据
data_path = 'GSE146249(TNBC)\data\GSE169246_converted_for_python.h5ad'
adata = scv.read(data_path)

# 4. 数据清洗 (手动 Numpy 版)
print(f"清洗前维度: {adata.shape}")
if 'spliced' in adata.layers:
    spliced_mat = adata.layers['spliced']
    if hasattr(spliced_mat, 'sum'):
        spliced_counts = spliced_mat.sum(axis=1)
        if hasattr(spliced_counts, 'A1'): 
            spliced_counts = spliced_counts.A1
        else:
            spliced_counts = np.array(spliced_counts).flatten()
    else:
        spliced_counts = np.sum(spliced_mat, axis=1)

    spliced_counts = np.nan_to_num(spliced_counts)
    cell_mask = spliced_counts > 0
    adata = adata[cell_mask].copy()
    print(f"清洗后维度: {adata.shape}")
else:
    print("警告: 未找到 'spliced' 层，跳过清洗。")

# 5. 预处理
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# 6. 配置参数 (根据您提供的 JSON 结构)
cfg = dv.Constants.default_configs.copy()

# CPU 适配修改
cfg['trainer']['epochs'] = 100  # CPU 跑慢，先设 30 轮
cfg['n_gpu'] = 0              # 显式关闭 GPU
cfg['loss']['args']['inner_batch_size'] = 64 # 降低内存压力
cfg['trainer']['save_dir'] = os.path.join(save_dir, 'checkpoints')

print("开始训练 DeepVelo (VeloGCN 模式)...")
# 7. 训练与自动预测
# 注意：dv.train 会自动完成训练、预测并将 'velocity' 层写入 adata
trainer = dv.train(adata, cfg)

# 10. 保存数据
out_path = 'GSE146249(TNBC)/scripts/deepvelo_velocity/MonoMacro_Part1_deepvelo.h5ad'
adata.write(out_path, compression='gzip')

# 8. 下游分析与绘图
# 直接使用 adata，因为 velocity 已经存在了
print("构建 Velocity Graph...")
scv.tl.velocity_graph(adata, n_jobs=1)

print("计算 Latent Time...")
scv.tl.velocity_pseudotime(adata)

# 标签处理
label_key = 'condition'
if label_key not in adata.obs.columns:
    sc.tl.leiden(adata)
    label_key = 'leiden'

# 9. 绘图保存
print("正在保存图片...")
# === 强制将 UMAP 坐标转为 Numpy 格式 ===
if isinstance(adata.obsm['UMAP'], pd.DataFrame):
    print("检测到 UMAP 为 DataFrame，正在转换为 Numpy array...")
    adata.obsm['UMAP'] = adata.obsm['UMAP'].to_numpy()
if 'PCA' in adata.obsm and isinstance(adata.obsm['PCA'], pd.DataFrame):
    adata.obsm['PCA'] = adata.obsm['PCA'].to_numpy()
# ======================================================

    
scv.pl.velocity_embedding_stream(
    adata, basis='UMAP', color=label_key,
    title='DeepVelo: Sensitive -> Resistant',
    show=False
)
plt.savefig(os.path.join(save_dir, 'deepvelo_stream.png'), dpi=300, bbox_inches='tight')
plt.close()

# Latent Time
scv.pl.scatter(
    adata, color='velocity_pseudotime', color_map='gnuplot', size=80,
    title='DeepVelo Latent Time', basis='UMAP',
    show=False
)
plt.savefig(os.path.join(save_dir, 'deepvelo_pseudo_time.png'), dpi=300, bbox_inches='tight')
plt.close()


print(f"DeepVelo 分析全部完成，结果保存至: {out_path}")