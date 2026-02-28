import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Set Chinese font (if needed)
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

# ==========================================
# 1. Core calculation functions
# ==========================================
def calculate_direction_scores(adata, source_group, target_group, 
                               group_key='condition', velocity_layer='velo_GemOut', expression_layer='M_total'):
    # Extract data
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
    # Calculate target centroid
    centroid = np.mean(target_X, axis=0)
    
    # Calculate cosine similarity
    cos_sims = []
    epsilon = 1e-8
    for i in range(source_X.shape[0]):
        velocity_vec = source_V[i]
        ideal_vec = centroid - source_X[i]  # Vector = end point - start point
        dot_product = np.dot(velocity_vec, ideal_vec)
        norm_v = np.linalg.norm(velocity_vec)
        norm_ideal = np.linalg.norm(ideal_vec)
        sim = dot_product / (norm_v * norm_ideal + epsilon) if norm_v > 0 else 0
        cos_sims.append(sim)
        
    return np.array(cos_sims)

# ==========================================
# 2. Calculate two sets of data and merge
# ==========================================
# Load data
adata = sc.read_h5ad("SigflowVelo/benchmark/GSE164551(MM)/data/SigflowVelo_data/Sen2Res.h5ad") 

# Positive signal: Sensitive -> Resistant
scores_signal = calculate_direction_scores(adata, 'sensitive', 'resistant')

# Negative control: Resistant -> Sensitive
scores_control = calculate_direction_scores(adata, 'resistant', 'sensitive')

# Merge all data points
all_scores = np.concatenate([scores_signal, scores_control])

# ==========================================
# 3. Draw single violin plot containing all data points (yellow and purple scatter points)
# ==========================================
# Convert to angles (0-180 degrees)
angles = np.degrees(np.arccos(all_scores))
print(f"Maximum angle: {np.max(angles)}")
print(f"Minimum angle: {np.min(angles)}")

# ==========================================
# 4. Draw rose plot (lower half plane)
# ==========================================
width_per_hist = 15
# Bin angles (every 15 degrees)
bins = np.arange(0, 181, width_per_hist)  # 0-180 degrees, every 15 degrees
hist, bin_edges = np.histogram(angles, bins=bins)

# Convert to percentages
total_cells = len(angles)
percentages = (hist / total_cells) * 100  # Percentages

# Create semi-circle rose plot
fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'projection': 'polar'})

# Calculate each bar position (center angle, converted to radians)
theta = np.deg2rad(bin_edges[:-1] + width_per_hist / 2)  # Center angle of each interval
width = np.deg2rad(width_per_hist)  # Bar width (radians)

# Draw bars (only 0-180 degrees)
bars = ax.bar(theta, percentages, width=width, alpha=0.8,
              edgecolor='white', linewidth=1.5)

# Set color mapping: smaller angles (better alignment) have warmer colors
bin_centers = bin_edges[:-1] + width_per_hist / 2  # Center angle of each interval (0-180)
# Normalize angles to [0,1] range, 0° corresponds to 1 (deepest), 180° corresponds to 0 (lightest)
color_norm = 1 - (bin_centers / 180.0)  # Invert so 0° has deepest color
colors = plt.cm.YlOrRd(color_norm)
for bar, color in zip(bars, colors):
    bar.set_facecolor(color)

# Set semi-circle display
ax.set_thetamin(0)
ax.set_thetamax(180)
ax.set_theta_zero_location('E')  # 0 degrees on the right
ax.set_theta_direction(-1)  # Clockwise direction

# Set radial grid and labels - fixed to 0-40%
ax.set_ylim(0, 50)  # Fixed radial range 0-40%
ax.set_rlabel_position(0)
ax.yaxis.grid(True, alpha=0.3)

# Set radial ticks (fixed at 0, 10, 20, 30, 40%)
ticks = [0, 10, 20, 30, 40, 50]
ax.set_yticks(ticks)
ax.set_yticklabels([f'{tick:.0f}%' for tick in ticks])

# Add angle labels
angle_labels = [f'{int(b)}°' for b in bin_edges]
ax.set_xticks(np.deg2rad(bin_edges))
ax.set_xticklabels(angle_labels)

# Add statistical information
median_angle = np.median(angles)
median_cosine = np.median(all_scores)

# Add statistical information in the center
ax.text(10, 7,
        f'Median angle = {median_angle:.1f}°\n'
        f'Median cosine = {median_cosine:.3f}\n',
        ha='right', va='bottom', fontsize=15,
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

plt.title('Rose Plot of Velocity Direction',
          fontsize=14, pad=20)

plt.tight_layout()
save_path = "output/figures/Roseplot_Cos_FlowVelo.svg"
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.show()

# Output save path
print(f"Rose plot saved to: {save_path}")