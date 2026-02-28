import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv
from scipy.stats import ks_2samp  # Import ks_2samp function
from evaluation_Metric import inner_cluster_coh
from scipy.sparse import csr_matrix
from evaluation_Metric import cross_boundary_correctness

# 1. Load data
adata_Sen2Res = sc.read_h5ad("SigflowVelo/benchmark/GSE164551(MM)/data/SigflowVelo_data/Sen2Res.h5ad")
print(adata_Sen2Res)
# Handle NaN in velo_hat
# adata_Sen2Res.layers['velo_hat'][np.isnan(adata_Sen2Res.layers['velo_hat'])] = 0
# print("NaN replaced with 0")

# Set parameters
velocity_pseudotime = 'latent_time' # Pseudotime column name
condition_key = 'condition'  # Condition column name - renamed to condition_key for clarity
sensitive = 'sensitive'  # Sensitive cell group name
resistant = 'resistant'  # Resistant cell group name

# Check and calculate KS statistic
print(f"\nCalculating KS Statistic (Latent Time Separation) ---")

if velocity_pseudotime not in adata_Sen2Res.obs.columns:
    print(f"Error: '{velocity_pseudotime}' not found in adata.obs.")
else:
    # Extract Source and Target times
    t_source = adata_Sen2Res.obs[adata_Sen2Res.obs[condition_key] == sensitive][velocity_pseudotime].values
    t_target = adata_Sen2Res.obs[adata_Sen2Res.obs[condition_key] == resistant][velocity_pseudotime].values
    
    if len(t_source) == 0 or len(t_target) == 0:
        print(f"Error: One of the groups ('{sensitive}' or '{resistant}') is empty in column '{condition_key}'.")
    else:
        # Calculate KS
        # ks_stat, p_val = ks_2samp(t_source, t_target)
        ks_statistic, p_value = ks_2samp(t_source, t_target, alternative="greater")
        # Check directionality (Target time should be later than Source)
        mean_source = np.mean(t_source)
        mean_target = np.mean(t_target)
        direction_correct = mean_target > mean_source
        
        print(f" -> KS Statistic: {ks_statistic:.4f} (Ideal ~ 1.0)")
        print(f" -> P-value: {p_value:.4e}")
        print(f" -> Mean Latent Time: {sensitive}={mean_source:.3f}, {resistant}={mean_target:.3f}")
        print(f" -> Time Direction Correct? {direction_correct}")