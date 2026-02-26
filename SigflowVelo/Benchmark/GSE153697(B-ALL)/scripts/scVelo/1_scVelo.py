import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')

plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
scv.set_figure_params(dpi=100, dpi_save=300, format='png', transparent=True)


def run_scvelo_analysis(adata_path, output_dir="scvelo_results",
                        min_shared_counts=20, min_shared_cells=10,
                        n_neighbors=30, n_pcs=30, n_top_genes=2000):
    import os
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 70)
    print("开始scVelo velocity分析")
    print("=" * 70)

    # Load data
    adata = sc.read_h5ad(adata_path)
    print(f"adata.shape: {adata.shape}")
    print(f"layers: {list(adata.layers.keys())}")

    # Check whether the data contains the required layers for velocity
    if 'spliced' not in adata.layers or 'unspliced' not in adata.layers:
        print("error: spliced or unspliced layer missing")
        return None

    # set obs and var
    if 'n_counts' not in adata.obs:
        adata.obs['n_counts'] = adata.layers['spliced'].sum(axis=1).A1 if hasattr(adata.layers['spliced'], 'A1') else \
        adata.layers['spliced'].sum(axis=1)

    # 4. Data preprocessing

    # Check if there are any high-mutation genes and PCA present
    if 'highly_variable' not in adata.var:
        print(" Data preprocessing...")
        sc.pp.filter_genes(adata, min_cells=min_shared_cells)
        sc.pp.filter_cells(adata, min_genes=200)

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='seurat')

    # 检查是否已经有PCA
    if 'X_pca' not in adata.obsm:
        print("  Calculating PCA...")
        scv.pp.pca(adata, n_comps=n_pcs)

    # Calculate moments
    scv.pp.moments(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    # Calculate RNA velocity
    try:
        print("  Trying dynamic model...")
        scv.tl.recover_dynamics(adata, n_jobs=4)
        scv.tl.velocity(adata, mode='dynamical')
        velocity_mode = 'dynamical'
        print("  ✓ dynamic model")
    except Exception as e:
        print(f"  Failed dymanic model: {e}")
        print("  Trying deterministic model...")
        scv.tl.velocity(adata, mode='deterministic')
        velocity_mode = 'deterministic'
        print("  ✓ deterministic model")

    # Calculating velocity graph and embedding
    print("\n7. Calculating velocity graph...")
    scv.tl.velocity_graph(adata)

    # 10. Calculating latent time（for dynamic model）
    if velocity_mode == 'dynamical':
        print("\n10. Calculating latent time...")
        scv.tl.latent_time(adata)

    # save adata
    output_adata_path = os.path.join(output_dir, 'adata_with_scVelo.h5ad')
    adata.write_h5ad(output_adata_path, compression='gzip')
    return adata


def main():
    print("RNA Velocity analysis")
    print("=" * 60)


    # Run scVelo
    adata_with_velocity = run_scvelo_analysis(
        adata_path="../../data/scPDAC_with_su_normalized.h5ad",
        output_dir="../../results/scVelo/adata",
        min_shared_counts=20,
        min_shared_cells=10,
        n_neighbors=30,
        n_pcs=30,
        n_top_genes=2000
    )


if __name__ == "__main__":
    main()