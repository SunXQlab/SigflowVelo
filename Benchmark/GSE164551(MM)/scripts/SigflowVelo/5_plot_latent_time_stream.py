import scanpy as sc
import scvelo as scv
import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
from scipy.stats import spearmanr
warnings.filterwarnings('ignore')


def _ensure_umap_numpy(adata, key='X_umap'):
    # scVelo plotting expects ndarray-like slicing for embeddings
    if key not in adata.obsm:
        raise KeyError(f"adata.obsm does not contain '{key}'. Cannot plot with basis='umap'.")
    if isinstance(adata.obsm[key], pd.DataFrame):
        adata.obsm[key] = adata.obsm[key].to_numpy()
    else:
        adata.obsm[key] = np.asarray(adata.obsm[key])


def main():
    # Read input data from previous step
    adata = sc.read("SigflowVelo/benchmark/GSE164551(MM)/data/SigflowVelo_data/Sen2Res.h5ad")
    print(adata)

    # Windows-friendly + reproducible
    scv.settings.n_jobs = 1

    # Ensure UMAP is numpy to avoid pandas slicing errors in scVelo plotting
    _ensure_umap_numpy(adata, key="X_umap")

    # ------------------------------------------------------------
    # Shared configuration (as in your scripts)
    # ------------------------------------------------------------
    vkey = "velo_GemOut"
    xkey = "Imputate"

    # ------------------------------------------------------------
    # Build neighbors once (required for velocity_graph in general)
    # (Keeps behavior consistent and avoids relying on precomputed neighbors.)
    # ------------------------------------------------------------
    sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_umap", random_state=0)

    # ------------------------------------------------------------
    # 1) Velocity graph (used by both downstream streamline 
    # ------------------------------------------------------------
    scv.tl.velocity_graph(adata, vkey=vkey, xkey=xkey)

    # ------------------------------------------------------------
    # 2) Streamline plot (output path adjusted)
    # ------------------------------------------------------------
    # Compute velocity embedding in UMAP space
    scv.tl.velocity_embedding(adata, basis="umap", vkey=vkey)

    ax_stream = scv.pl.velocity_embedding_stream(
    adata,
    color="condition",
    density=2,
    smooth=0.5,
    title="RNA Velocity Stream Plot Of SigflowVelo",
    show=False,
    vkey=vkey,
    linewidth=0.5,
    arrow_size=0.5,
    legend_loc="right",
    legend_fontsize=8,
    )

    fig_stream = ax_stream.figure if hasattr(ax_stream, "figure") else plt.gcf()
    fig_stream.canvas.draw()
    fig_stream.savefig(
        "output/figures/streamline_prediction.svg",
        bbox_inches="tight",
        pad_inches=0.1,
        transparent=True,
        dpi=300,
    )
    plt.close(fig_stream)

    # ------------------------------------------------------------
    # 3) Latent time dot plot (output path adjusted)
    #    Note: this assumes 'latent_time' already exists in adata.obs (as in your script).
    # ------------------------------------------------------------
    if "latent_time" not in adata.obs:
        raise KeyError(
            "adata.obs does not contain 'latent_time'. "
            "Your original script plotted an existing latent_time; "
            "if you want to compute latent_time here, tell me your dynamical setup."
        )

    ax_lt = scv.pl.scatter(
        adata,
        color="latent_time",
        color_map="plasma",
        size=200,
        title="Latent Time on UMAP ",
        show=False,
    )
    fig_lt = ax_lt.figure if hasattr(ax_lt, "figure") else plt.gcf()
    fig_lt.canvas.draw()
    fig_lt.savefig(
        "output/figures/umap_of_latent_time.svg",
        format="svg",
        transparent=True,
        dpi=300,
    )
    plt.close(fig_lt)

    # ------------------------------------------------------------
    # 4) Save outputs (save to output/data/ directory)
    # ------------------------------------------------------------
    adata.write_h5ad("SigflowVelo/benchmark/GSE164551(MM)/data/SigflowVelo_data/adata_with_Flowsig.h5ad")
    adata.write("SigflowVelo/benchmark/GSE164551(MM)/data/SigflowVelo_data/Sen2Res_with_velocity_pseudotime.h5ad")
    print("Spearman r =", spearmanr(adata.obs["latent_time"], adata.obs["velo_GemOut_pseudotime"], nan_policy="omit").statistic)


if __name__ == "__main__":
    mp.freeze_support()
    main()