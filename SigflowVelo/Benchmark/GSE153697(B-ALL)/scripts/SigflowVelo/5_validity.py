import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics.pairwise import cosine_similarity
from scipy import sparse

def sensitive_cosines_calculated(adata, group_key='group', velocity_key='velo_GemOut_umap'):
    velocity_2d = adata.obsm[velocity_key]  # 2D velocity vectors
    groups = adata.obs[group_key]  # Cell groups
    umap_coords = adata.obsm['X_umap']  # UMAP coordinates

    # Calculate center positions of sensitive and resistant cells
    sensitive_mask = groups == 'sensitive'
    resistant_mask = groups == 'resistant'

    sensitive_center = umap_coords[sensitive_mask].mean(axis=0)
    resistant_center = umap_coords[resistant_mask].mean(axis=0)

    # For each sensitive cell: calculate the vector pointing from this cell to the center of the resistant cell
    if np.sum(sensitive_mask) > 0:
        sensitive_cells_coords = umap_coords[sensitive_mask]  # Coordinates of all sensitive cells
        sensitive_velocity = velocity_2d[sensitive_mask]  # 2D velocity vectors of all sensitive cells

        # Calculate the direction vector each sensitive cell points towards the center of the resistant cells
        sensitive_to_resistant_vectors = resistant_center - sensitive_cells_coords

        # Calculate magnitude of each direction vector
        direction_norms = np.linalg.norm(sensitive_to_resistant_vectors, axis=1)
        # Calculate magnitude of each velocity vector
        velocity_norms = np.linalg.norm(sensitive_velocity, axis=1)

        # Calculate dot product: velocity vector of each sensitive cell with corresponding direction vector
        dot_products = np.sum(sensitive_velocity * sensitive_to_resistant_vectors, axis=1)

        # Calculate cosine value: dot product / (direction magnitude * velocity magnitude)
        sensitive_cosines = np.zeros(len(sensitive_cells_coords))
        valid_mask = (direction_norms > 0) & (velocity_norms > 0)
        sensitive_cosines[valid_mask] = dot_products[valid_mask] / (
                    direction_norms[valid_mask] * velocity_norms[valid_mask])
        return (sensitive_cosines)

def resistant_cosines_calculated(adata, group_key='group', velocity_key='velo_GemOut_umap'):
    velocity_2d = adata.obsm[velocity_key]  # 2D velocity vectors
    groups = adata.obs[group_key]  # Cell groups
    umap_coords = adata.obsm['X_umap']  # UMAP coordinates

    # Calculate center positions of sensitive and resistant cells
    sensitive_mask = groups == 'sensitive'
    resistant_mask = groups == 'resistant'

    sensitive_center = umap_coords[sensitive_mask].mean(axis=0)
    resistant_center = umap_coords[resistant_mask].mean(axis=0)

    # For each resistant cell: calculate the vector pointing from the center of the resistant cell to this cell
    if np.sum(resistant_mask) > 0:
        resistant_cells_coords = umap_coords[resistant_mask]  # Coordinates of all resistant cells
        resistant_velocity = velocity_2d[resistant_mask]  # 2D velocity vectors of all resistant cells

        # Calculate the direction vector the center of the sensitive cells points towards each resistant cell
        sensitive_to_resistant_vectors = resistant_cells_coords - sensitive_center

        # Calculate magnitude of each direction vector
        direction_norms = np.linalg.norm(sensitive_to_resistant_vectors, axis=1)
        # Calculate magnitude of each velocity vector
        velocity_norms = np.linalg.norm(resistant_velocity, axis=1)

        # Calculate dot product: velocity vector of each resistant cell with corresponding direction vector
        dot_products = np.sum(resistant_velocity * sensitive_to_resistant_vectors, axis=1)

        # Calculate cosine value: dot product / (direction magnitude * velocity magnitude)
        resistant_cosines = np.zeros(len(resistant_cells_coords))
        valid_mask = (direction_norms > 0) & (velocity_norms > 0)
        resistant_cosines[valid_mask] = dot_products[valid_mask] / (
                direction_norms[valid_mask] * velocity_norms[valid_mask])
        return (resistant_cosines)

def rose_plot_of_velocity_direction(adata, group_key, velocity_key):
    sensitive_cosines = sensitive_cosines_calculated(adata, group_key, velocity_key)
    resistant_cosines = resistant_cosines_calculated(adata, group_key, velocity_key)
    score_total = np.concatenate((sensitive_cosines, resistant_cosines))
    # Convert to angles (0-180 degrees)
    angles = np.degrees(np.arccos(score_total))
    width_per_hist = 15
    # Bin angles (every 15 degrees as an interval)
    bins = np.arange(0, 181, width_per_hist)  # 0-180degrees，15 degrees per bin
    hist, bin_edges = np.histogram(angles, bins=bins)

    # Convert to percentages
    total_cells = len(angles)
    print(total_cells)
    percentages = (hist / total_cells) * 100

    # Create semicircular rose plot
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'projection': 'polar'})

    # Calculate position of each bar (center angle, convert to radians)
    theta = np.deg2rad(bin_edges[:-1] + width_per_hist / 2)  # central angle for every bin
    width = np.deg2rad(width_per_hist)

    # Draw bars
    bars = ax.bar(theta, percentages, width=width, alpha=0.8,
                  edgecolor='white', linewidth=1.5)

    # Set color mapping: smaller angles (better alignment) have warmer colors
    bin_centers = bin_edges[:-1] + width_per_hist / 2  # Center angle of every bin（0-180）
    color_norm = 1 - (bin_centers / 180.0)
    colors = plt.cm.YlOrRd(color_norm)
    for bar, color in zip(bars, colors):
        bar.set_facecolor(color)

    # Set semicircular display
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    ax.set_theta_zero_location('E')
    ax.set_theta_direction(-1)

    # Set radial grid
    ax.set_rlabel_position(0)
    ax.yaxis.grid(True, alpha=0.3)
    max_percentage = max(percentages)
    tick_step = 10
    max_tick = int(np.ceil(max_percentage / tick_step) * tick_step)
    ticks = np.arange(0, max_tick + 1, tick_step)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{tick:.0f}%' for tick in ticks])

    # Add angle labels
    angle_labels = [f'{int(b)}°' for b in bin_edges]
    ax.set_xticks(np.deg2rad(bin_edges))
    ax.set_xticklabels(angle_labels)

    plt.title('rose plot of velocity direction',
              fontsize=14, pad=20)
    plt.tight_layout()
    save_path = "../../results/SigflowVelo/plot/roseplot_of_cosine.svg"
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

#KS statistic for latent time of two group
def Kolmogorov_Smirnov_divergence(adata, group_key='group', latent_time_key='latent_time'):
    sensitive_mask = adata.obs[group_key] == 'sensitive'
    resistant_mask = adata.obs[group_key] == 'resistant'
    sensitive_times = adata.obs.loc[sensitive_mask, latent_time_key].to_numpy().astype(float)
    resistant_times = adata.obs.loc[resistant_mask, latent_time_key].to_numpy().astype(float)
    # Calculate KS statistic
    ks_statistic, p_value = stats.ks_2samp(sensitive_times, resistant_times, alternative="greater")
    print(f"ks_statistic:{ks_statistic}")
    print (f"p_value:{p_value}")

if __name__ == "__main__":
    adata = sc.read_h5ad("../../results/SigflowVelo/adata/adata_with_embedded_velocity.h5ad")
    rose_plot_of_velocity_direction(adata, group_key='group', velocity_key='velo_GemOut_umap')
    Kolmogorov_Smirnov_divergence(adata, group_key='group', latent_time_key='latent_time')
