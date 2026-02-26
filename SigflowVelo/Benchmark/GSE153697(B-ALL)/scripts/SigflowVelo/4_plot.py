import torch
import scanpy as sc
import pandas as pd
# import TFvelo as TFv
import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
import warnings
import time
from tqdm import tqdm
from collections import Counter
import scvelo as scv
import seaborn as sns
from scipy.stats import spearmanr

def plot_for_streamline_colored_by_group(adata_copy):
    scv.settings.verbosity = 0  # 关闭冗余警告
    scv.tl.velocity_graph(adata_copy, basis='X_umap', vkey='velo_GemOut', xkey='spliced')
    scv.pl.velocity_embedding_stream(
        adata_copy,
        basis='X_umap',
        vkey='velo_GemOut',
        color='group',
        density=2,
        smooth=0.5,
        linewidth=0.5,
        arrow_size=0.5,
        legend_loc='right',
        legend_fontsize=8,
        title='Inferred RNA Velocity Streamline',
        show=False
    )

    plt.savefig("../../resluts/SigflowVelo/plot/streamline_prediction_colored_by_group.svg",
                format='svg',
                bbox_inches='tight',
                pad_inches=0.1,
                transparent=True,
                dpi=300)
    plt.show()
    plt.close()

def umapplot_by_latent_time(adata_copy):
    plt.figure(figsize=(7, 5))
    scv.pl.scatter(
        adata_copy,
        color='latent_time',
        color_map='gnuplot',
        size=20,
        title='Latent Time',
        show=False
    )
    plt.savefig("../../resluts/SigflowVelo/plot/umap_by_latent_time.svg",
                format='svg',
                bbox_inches='tight',
                pad_inches=0.1,
                transparent=True,
                dpi=300)
    plt.show()
    plt.close()
    scv.pl.scatter(
        adata_copy,
        color='velo_GemOut_pseudotime',
        color_map='gnuplot',
        size=20,
        title='pseudotime',
        show=False
    )
    plt.show()

    correlation = spearmanr(adata_copy.obs['latent_time'], adata_copy.obs['velo_GemOut_pseudotime'])[0]
    print(f"spearman_correlation: {correlation}")

if __name__ == '__main__':
    adata_Sen2Res = sc.read_h5ad("../../resluts/SigflowVelo/adata/adata_with_calculated_velocity.h5ad")
    adata_Sen2Res.layers['spliced'] = adata_Sen2Res.X.copy()
    print (adata_Sen2Res)
    umapplot_by_latent_time(adata_Sen2Res)
    plot_for_streamline_colored_by_group(adata_Sen2Res)
    adata_Sen2Res.write_h5ad("../../resluts/SigflowVelo/adata/adata_with_embedded_velocity.h5ad")