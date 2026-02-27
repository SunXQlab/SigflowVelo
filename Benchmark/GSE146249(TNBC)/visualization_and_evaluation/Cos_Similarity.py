import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# 设置中文字体（如果需要）
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

# ==========================================
# 1. 核心计算函数
# ==========================================
def calculate_direction_scores(adata, source_group, target_group, 
                               group_key='condition', velocity_layer='velo_GemOut', expression_layer='M_total'):
    # 提取数据
    if expression_layer in adata.layers:
        X = adata.layers[expression_layer]
    else:
        X = adata.X
    V = adata.layers[velocity_layer]
    
    source_mask = adata.obs[group_key] == source_group
    target_mask = adata.obs[group_key] == target_group
    
    source_X = X[source_mask]
    source_V = V[source_mask]
    target_X = X[target_mask]
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
# 2. 计算两组数据并合并
# ==========================================
# 加载数据
adata = sc.read_h5ad("script_5_evaluation\Mono-Macro_flowsig_velocity_results.h5ad") 

# 正向信号: Sensitive -> Resistant
scores_signal = calculate_direction_scores(adata, 'sensitive', 'resistant')

# 负向对照: Resistant -> Sensitive
scores_control = calculate_direction_scores(adata, 'resistant', 'sensitive')

# 合并所有数据点
all_scores = np.concatenate([scores_signal, scores_control])

# ==========================================
# 3. 绘制包含所有数据点的单个小提琴图（黄色和紫色散点）
# ==========================================
plt.figure(figsize=(8, 10))

# 颜色设置
main_color = '#2E8B57'  # 小提琴图颜色（海绿色）
signal_color = '#FFD700'  # 黄色，代表正向信号
control_color = '#800080'  # 紫色，代表反向对照

# 准备数据（单个位置）
data = [all_scores]
position = [1]
label = 'All Directions'

def rose_plot_of_velocity_direction(scores_signal, scores_control):
    sensitive_cosines = scores_signal
    resistant_cosines = scores_control
    score_total = np.concatenate((sensitive_cosines, resistant_cosines))
    # 转换为角度（0-180度）
    angles = np.degrees(np.arccos(score_total))

    width_per_hist=9
    # 将角度分箱（每9度一个区间）
    bins = np.arange(0, 181, width_per_hist)  # 0-180度，每9度一个区间
    hist, bin_edges = np.histogram(angles, bins=bins)

    # 创建半圆玫瑰图
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'projection': 'polar'})

    # 计算每个条形的位置（中心角度，转换为弧度）
    theta = np.deg2rad(bin_edges[:-1] + width_per_hist/2)  # 每个区间的中心角度
    width = np.deg2rad(width_per_hist)  # 条形的宽度（弧度）

    # 绘制条形（只绘制0-180度）
    bars = ax.bar(theta, hist, width=width, alpha=0.8,
                  edgecolor='white', linewidth=1.5)

    # 设置颜色映射：角度越小（对齐越好）颜色越暖
    bin_centers = bin_edges[:-1] +width_per_hist/2 # 每个区间的中心角度（0-180）
    # 归一化角度到[0,1]范围，0°对应1（最深），180°对应0（最浅）
    color_norm = 1 - (bin_centers / 180.0)  # 反转，使0°颜色最深
    colors = plt.cm.YlOrRd(color_norm)
    for bar, color in zip(bars, colors):
        bar.set_facecolor(color)

    # 设置半圆显示
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    ax.set_theta_zero_location('E')  # 0度在右侧
    ax.set_theta_direction(-1)  # 顺时针方向

    # 设置径向网格
    ax.set_rlabel_position(0)  # 径向标签位置
    ax.yaxis.grid(True, alpha=0.3)

    # 添加角度标签
    angle_labels = [f'{int(b)}°' for b in bin_edges]
    ax.set_xticks(np.deg2rad(bin_edges))
    ax.set_xticklabels(angle_labels)

    # 添加统计信息
    median_angle = np.median(angles)
    median_cosine = np.median(score_total)

    # 在中心添加统计信息
    ax.text(1000, 700,
            f'Median angle = {median_angle:.1f}°\n'
            f'Median cosine = {median_cosine:.3f}\n',
            ha='right', va='bottom', fontsize=15,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    plt.title('rose plot of velocity direction',
              fontsize=14, pad=20)

    plt.tight_layout()
    plt.show()

rose_plot_of_velocity_direction(scores_signal, scores_control)