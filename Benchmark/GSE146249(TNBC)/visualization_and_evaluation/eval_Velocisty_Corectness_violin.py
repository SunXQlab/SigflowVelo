import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import sparse

# 设置中文字体（保留原设置）
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

# ==========================================
# 0. 配置区域 (请在此处修改文件路径和图层名称)
# ==========================================
model_configs = {
    "FlowVelo": {
        "path": "Mono-Macro_FSvelo_results.h5ad",  # FSvelo的数据路径
        "velocity_layer": "velo_GemOut",    # 请修改为FSvelo对应的layer名，如 'velocity_FS'
        "expression_layer": "Imputate"       # 表达量图层
    },
    "TFVelo": {
        "path": "Mono-Macro_TFvelo_result.h5ad",  # TFvelo的数据路径
        "velocity_layer": "velocity",        # 请修改为TFvelo对应的layer名
        "expression_layer": "M_total"
    }
}

# 定义组别信息
GROUP_KEY = 'condition'
SOURCE_GROUP = 'sensitive'
TARGET_GROUP = 'resistant'

# ==========================================
# 1. 核心计算函数 (保持不变)
# ==========================================
def calculate_direction_scores(adata, source_group, target_group, 
                               group_key='condition', velocity_layer='velo_GemOut', expression_layer='M_total'):
    # 提取数据
    if expression_layer in adata.layers:
        X = adata.layers[expression_layer]
    else:
        X = adata.X
        
    # 检查velocity layer是否存在
    if velocity_layer not in adata.layers:
        raise ValueError(f"图层 '{velocity_layer}' 不在 adata.layers 中，请检查配置。")
        
    V = adata.layers[velocity_layer]
    
    source_mask = adata.obs[group_key] == source_group
    target_mask = adata.obs[group_key] == target_group
    
    source_X = X[source_mask]
    source_V = V[source_mask]
    target_X = X[target_mask]
# --- 1. 稀疏矩阵处理 ---
    # 为了后续方便检测NaN和计算，如果是稀疏矩阵，先转为稠密矩阵
    if sparse.issparse(source_X): source_X = source_X.toarray()
    if sparse.issparse(source_V): source_V = source_V.toarray()
    if sparse.issparse(target_X): target_X = target_X.toarray()

    # --- 2. NaN 值预处理 (按列/基因剔除) ---
    # 检查 source_V, source_X, target_X 中哪些列(基因)包含 NaN
    # 只要某个基因在任意一个细胞中是 NaN，为了保证维度一致和计算准确，我们都在本次计算中剔除该基因
    
    nan_mask_V = np.isnan(source_V).any(axis=0)
    nan_mask_X_src = np.isnan(source_X).any(axis=0)
    nan_mask_X_tgt = np.isnan(target_X).any(axis=0)
    
    # 合并掩码：任何一个矩阵中出问题的基因都标记为 True
    genes_to_drop_mask = nan_mask_V | nan_mask_X_src | nan_mask_X_tgt
    
    n_dropped = np.sum(genes_to_drop_mask)
    if n_dropped > 0:
        print(f"  [NaN处理] 检测到 {n_dropped} 个基因包含 NaN 值，已从计算中剔除这些列。")
        # 保留正常的基因
        valid_genes = ~genes_to_drop_mask
        source_X = source_X[:, valid_genes]
        source_V = source_V[:, valid_genes]
        target_X = target_X[:, valid_genes]
    
    # 检查是否还有剩余基因
    if source_X.shape[1] == 0:
        print("  [错误] 剔除 NaN 后没有剩余的有效基因，无法计算。")
        return np.array([])
        
    # 检查是否有剩余细胞 (通常不会变，除非本身为空)
    if source_X.shape[0] == 0 or target_X.shape[0] == 0:
        print("  [错误] 有效细胞数量为 0。")
        return np.array([])
    
    # 计算目标质心
    centroid = np.mean(target_X, axis=0)
    
    # 计算余弦相似度
    cos_sims = []
    epsilon = 1e-8
    for i in range(source_X.shape[0]):
        velocity_vec = source_V[i]
        ideal_vec = centroid - source_X[i] # 向量 = 终点 - 起点
        dot_product = np.dot(velocity_vec, ideal_vec)
        norm_v = np.linalg.norm(velocity_vec)
        norm_ideal = np.linalg.norm(ideal_vec)
        sim = dot_product / (norm_v * norm_ideal + epsilon) if norm_v > 0 else 0
        cos_sims.append(sim)
        
    return np.array(cos_sims)

# ==========================================
# 2. 循环加载数据并计算
# ==========================================
plot_data = [] # 存储用于绘图的数据结构
print("开始计算模型得分...")

for model_name, config in model_configs.items():
    print(f"正在处理模型: {model_name} ...")
    
    # 加载数据 (如果两个模型用同一个文件，实际应用中可以优化只读一次，这里为了灵活保持分别读取)
    adata = sc.read_h5ad(config["path"]) 
    
    # 正向信号: Sensitive -> Resistant
    scores_signal = calculate_direction_scores(
        adata, SOURCE_GROUP, TARGET_GROUP, 
        group_key=GROUP_KEY, 
        velocity_layer=config["velocity_layer"], 
        expression_layer=config["expression_layer"]
    )

    # 负向对照: Resistant -> Sensitive
    scores_control = calculate_direction_scores(
        adata, TARGET_GROUP, SOURCE_GROUP, 
        group_key=GROUP_KEY, 
        velocity_layer=config["velocity_layer"], 
        expression_layer=config["expression_layer"]
    )
    
    # 合并数据用于小提琴图整体轮廓
    all_scores = np.concatenate([scores_signal, scores_control])
    
    # 存储结果
    plot_data.append({
        "name": model_name,
        "all_scores": all_scores,
        "scores_signal": scores_signal,
        "scores_control": scores_control
    })

print("计算完成。")

# ==========================================
# 3. 绘制对比小提琴图
# ==========================================
plt.figure(figsize=(10, 8)) # 稍微加宽一点以容纳两个小提琴

# 颜色设置
main_color = '#2E8B57'   # 小提琴图轮廓颜色
signal_color = '#FFD700' # 黄色，正向
control_color = '#800080' # 紫色，反向

# 提取用于violinplot的数据列表和位置
violin_vectors = [d["all_scores"] for d in plot_data]
positions = [1, 2]
labels = [d["name"] for d in plot_data]

# 创建小提琴图
violin_parts = plt.violinplot(violin_vectors, positions=positions, 
                              showmeans=True, showmedians=True, 
                              showextrema=True,
                              widths=0.4,
                              points=1000)

# 自定义小提琴图样式
for pc in violin_parts['bodies']:
    pc.set_facecolor(main_color)
    pc.set_edgecolor(main_color)
    pc.set_alpha(0.8) # 稍微调低透明度，让散点更清楚
    pc.set_zorder(1)

# 设置统计线样式
violin_parts['cmeans'].set_color('blue')
violin_parts['cmeans'].set_linewidth(2)
violin_parts['cmedians'].set_color('black')
violin_parts['cmedians'].set_linewidth(2)
violin_parts['cmaxes'].set_color('gray')
violin_parts['cmaxes'].set_linewidth(1)
violin_parts['cmins'].set_color('gray')
violin_parts['cmins'].set_linewidth(1)

# 添加散点图（抖动点）
for i, d in enumerate(plot_data):
    pos = positions[i]
    scores_sig = d["scores_signal"]
    scores_ctr = d["scores_control"]
    
    # 抖动处理
    # 生成两组抖动，合并
    jitter_sig = np.random.normal(pos, 0.04, size=len(scores_sig))
    jitter_ctr = np.random.normal(pos, 0.04, size=len(scores_ctr))
    
    # 绘制正向信号点 (黄色)
    plt.scatter(jitter_sig, scores_sig, s=3, alpha=0.3, c=signal_color, 
                edgecolors='none', linewidth=0.3, zorder=0)
    
    # 绘制负向对照点 (紫色)
    plt.scatter(jitter_ctr, scores_ctr, s=3, alpha=0.3, c=control_color, 
                edgecolors='none', linewidth=0.3, zorder=0)

# 添加参考线
plt.axhline(y=0, color='gray', linestyle='--', linewidth=1.5, 
            alpha=0.7)

# 设置图形属性
plt.xticks(positions, labels, fontsize=13, fontweight='bold')
plt.ylabel('Velocity Correctness (Cosine Similarity)', fontsize=13, fontweight='bold')
plt.title('Velocity Directionality Comparison', fontsize=16, fontweight='bold', pad=15)
plt.ylim(-1.1, 1.1)
plt.grid(True, alpha=0.3, linestyle='--', axis='y')

# 添加图例 (只需添加一次，通用于两个模型)
legend_elements = [
    Line2D([0], [0], color=main_color, linewidth=4, alpha=0.5, label='Overall Distribution'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=signal_color, 
           markersize=10, label=f'{SOURCE_GROUP} → {TARGET_GROUP}'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=control_color, 
           markersize=10, label=f'{TARGET_GROUP} → {SOURCE_GROUP}'),
    Line2D([0], [0], color='blue', linestyle='-', linewidth=1.5, 
           label='Mean'),
    Line2D([0], [0], color='black', linestyle='-', linewidth=1.5, 
           label='Median'),
    Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, 
           label='Orthogonal (0.0)'),                      
]
plt.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.9)

# 调整布局
plt.tight_layout()

# 保存
save_path = "figures/Comparition_Velocity_Corectness_Violin.png"
# 确保目录存在 (可选)
import os
if not os.path.exists('figures'):
    os.makedirs('figures')
    
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.show()

print(f"对比图已保存至: {save_path}")
for d in plot_data:
    print(f"--- 模型: {d['name']} ---")
    print(f"  总细胞数: {len(d['all_scores'])}")
    print(f"  正向信号均值: {np.mean(d['scores_signal']):.4f}")
    print(f"  反向对照均值: {np.mean(d['scores_control']):.4f}")