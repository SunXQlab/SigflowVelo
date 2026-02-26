import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches

# 1. 准备数据 (请将此部分替换为您的真实数据读取代码，如 pd.read_csv('results.csv'))
data = {
    'Dataset': ['NK','NK','NK','NK','NK','Mono-Macro', 'Mono-Macro','Mono-Macro', 'Mono-Macro','Mono-Macro',
                'malignant','malignant','malignant','malignant','malignant'],
    'Model': ['FlowVelo', 'TFVelo','UniTVelo','DeepVelo','scVelo',
              'FlowVelo', 'TFVelo','UniTVelo','DeepVelo','scVelo',
              'FlowVelo','TFVelo','UniTVelo','DeepVelo','scVelo'],
    # 示例 KS 统计量数据
    'KS_statistic': [0.920,0.068, np.nan,np.nan,np.nan,0.495, 0.219,np.nan,np.nan,np.nan,
                     0.611,0.005,0.836,0.285,0],
    # 示例 P-value 数据
    'P_value': [5.3817e-89, 4.2475e-01,np.nan,np.nan,np.nan,9.94e-245,2.87e-46,np.nan,np.nan,np.nan,
                1.188e-144,0.977,5.030e-270,5.081e-32,1.0]
}
df = pd.DataFrame(data)
df['NegLog10_P'] = -np.log10(df['P_value'] + 1e-300)
# 2. 设置绘图风格
sns.set_theme(style="whitegrid")
plt.figure(figsize=(10, 6)) # 设置画布大小
fig, ax = plt.subplots(figsize=(10, 6))
# 3. 绘制气泡图 (Bubble Plot)
# x, y: 坐标轴
# size: 控制气泡大小的列名 (KS统计量)
# hue: 控制颜色的列名 (P值)
# sizes: (最小气泡尺寸, 最大气泡尺寸) -> 可根据需要调整
# palette: 颜色映射方案 (如 'viridis', 'Reds', 'Blues_r')
scatter = sns.scatterplot(
    data=df, 
    x='Dataset', 
    y='Model', 
    size='KS_statistic', 
    hue='NegLog10_P',
    sizes=(200, 1000), 
    palette='Reds',  # '_r' 表示颜色反转，通常P值越小越显著，看您希望深色代表大值还是小值
    alpha=0.8,            # 透明度
    edgecolor='black',    # 气泡边缘颜色
    linewidth=1.5
)

# 4. 美化图表
# ==========================================
# 核心修改：双图例分离法
# ==========================================

# 1. 获取所有句柄和标签
handles, labels = ax.get_legend_handles_labels()

# 2. 找到分割点 (KS_statistic 标题的位置)
try:
    split_idx = labels.index('KS_statistic')
except ValueError:
    split_idx = len(labels) // 2 

# 3. 拆分列表
# 注意：Seaborn生成的列表里，索引0通常是标题句柄，我们可以选择跳过它(切片[1:])自己加标题，或者保留
h_color = handles[:split_idx]
l_color = labels[:split_idx]
h_size = handles[split_idx:]
l_size = labels[split_idx:]

# 4. 创建第一个图例 (颜色)
# bbox_to_anchor=(1.02, 1) 表示放在图表右边缘紧贴着
legend_color = plt.legend(
    h_color[1:], l_color[1:],    # [1:] 去掉Seaborn自带的标题项，改用下面的 title 参数
    bbox_to_anchor=(1.02, 1), 
    loc='upper left', 
    borderaxespad=0.,
    title="-log10(P-value)",     # 自定义清晰的标题
    labelspacing=2,
    frameon=False                # 可选：去掉边框，让两列看起来更融合
)

# 【至关重要】将第一个图例手动“冻结”在图上
# 如果不加这行，第二次调用 plt.legend() 会自动清除掉前一个图例
ax.add_artist(legend_color)

# 5. 创建第二个图例 (大小)
# bbox_to_anchor=(1.25, 1) 表示放在第一个图例的更右边
# 您可以通过调整 1.25 这个数字来控制两列之间的距离
legend_size = plt.legend(
    h_size[1:], l_size[1:],      # [1:] 去掉自带标题
    bbox_to_anchor=(1.25, 1),    # 向右平移，形成第二列
    loc='upper left', 
    borderaxespad=0.,
    title="KS Statistic",
    labelspacing=2,
    frameon=False                # 去掉边框
)

plt.title('Model Evaluation: KS ststistic and P-value', fontsize=16, pad=20)
plt.xlabel('Dataset', fontsize=12)
plt.ylabel('Model', fontsize=12)

# 增加边距防止气泡被截断
plt.margins(0.2)
plt.tight_layout()

# 5. 保存或显示
plt.savefig('figures/Comparison_KS&p_bubble_tian&cao.png', dpi=300, bbox_inches='tight')
plt.show()