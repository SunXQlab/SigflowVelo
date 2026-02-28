import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os

plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False


def calculate_direction_scores(
    adata, source_group, target_group,
    group_key='condition', velocity_layer='velo_hat', expression_layer='M_total'
):
    # Expression
    X = adata.layers[expression_layer] if expression_layer in adata.layers else adata.X
    V = adata.layers[velocity_layer]

    source_mask = (adata.obs[group_key] == source_group).values if hasattr(adata.obs[group_key], "values") else (adata.obs[group_key] == source_group)
    target_mask = (adata.obs[group_key] == target_group).values if hasattr(adata.obs[group_key], "values") else (adata.obs[group_key] == target_group)

    source_X = X[source_mask]
    source_V = V[source_mask]
    target_X = X[target_mask]

    centroid = np.mean(target_X, axis=0)

    cos_sims = []
    eps = 1e-8
    for i in range(source_X.shape[0]):
        velocity_vec = source_V[i]
        ideal_vec = centroid - source_X[i]

        dot = np.dot(velocity_vec, ideal_vec)
        norm_v = np.linalg.norm(velocity_vec)
        norm_ideal = np.linalg.norm(ideal_vec)

        # When norm_v==0, set to 0 according to original logic (corresponding to 90°)
        sim = dot / (norm_v * norm_ideal + eps) if norm_v > 0 else 0.0
        cos_sims.append(sim)

    return np.asarray(cos_sims, dtype=np.float64)



# 1) Read data 
# Create output directory if it doesn't exist
os.makedirs("output/figures", exist_ok=True)

adata = sc.read_h5ad("SigflowVelo/benchmark/GSE164551(MM)/data/TFVelo_data/TFvelo.h5ad")

V = np.asarray(adata.layers["velo_hat"], dtype=np.float64)
nan_before = np.isnan(V).sum()
V = np.nan_to_num(V, nan=0.0, posinf=0.0, neginf=0.0)  # Replace only NaN/Inf (set Inf to 0 here)
adata.layers["velo_hat"] = V
print("NaN replaced in velo_hat:", int(nan_before), "->", int(np.isnan(adata.layers["velo_hat"]).sum()))


# 2) Calculate cosine similarity and convert to angles
scores_signal = calculate_direction_scores(adata, 'sensitive', 'resistant')
scores_control = calculate_direction_scores(adata, 'resistant', 'sensitive')
all_scores = np.concatenate([scores_signal, scores_control])

# Prevent numerical errors causing arccos input to exceed [-1,1]
all_scores = np.clip(all_scores, -1.0, 1.0)

angles = np.degrees(np.arccos(all_scores))
print("Maximum angle:", np.max(angles))
print("Minimum angle:", np.min(angles))


# 3) Official rose plot drawing (keep as original)
width_per_hist = 15
bins = np.arange(0, 181, width_per_hist)
hist, bin_edges = np.histogram(angles, bins=bins)

total_cells = len(angles)
percentages = (hist / total_cells) * 100

fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'projection': 'polar'})

theta = np.deg2rad(bin_edges[:-1] + width_per_hist/2)
width = np.deg2rad(width_per_hist)

bars = ax.bar(theta, percentages, width=width, alpha=0.8,
              edgecolor='white', linewidth=1.5)

bin_centers = bin_edges[:-1] + width_per_hist/2
color_norm = 1 - (bin_centers / 180.0)
colors = plt.cm.YlOrRd(color_norm)
for bar, color in zip(bars, colors):
    bar.set_facecolor(color)

ax.set_thetamin(0)
ax.set_thetamax(180)
ax.set_theta_zero_location('E')
ax.set_theta_direction(-1)

ax.set_ylim(0, 80)
ax.set_rlabel_position(0)
ax.yaxis.grid(True, alpha=0.3)

ticks = [0, 10, 20, 30, 40, 50, 60, 70, 80]
ax.set_yticks(ticks)
ax.set_yticklabels([f'{tick:.0f}%' for tick in ticks])

angle_labels = [f'{int(b)}°' for b in bin_edges]
ax.set_xticks(np.deg2rad(bin_edges))
ax.set_xticklabels(angle_labels)

median_angle = np.median(angles)
median_cosine = np.median(all_scores)

ax.text(10, 7,
        f'Median angle = {median_angle:.1f}°\n'
        f'Median cosine = {median_cosine:.3f}\n',
        ha='right', va='bottom', fontsize=15,
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

plt.title('Rose plot of velocity direction', fontsize=14, pad=20)

plt.tight_layout()
save_path = "output/figures/Roseplot_Cos_TFVelo.svg"
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.show()

print("Rose plot saved to:", save_path)