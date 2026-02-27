import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore')

if __name__ == "__main__":
    adata = sc.read("Mono-Macro_flowsig_velocity_results.h5ad")
    print (adata)
    # 计算velocity embedding（在UMAP空间中的velocity）
    print("  计算velocity embedding...")
    scv.tl.velocity_embedding(adata, basis='UMAP',vkey='velo_GemOut')
    scv.pl.velocity_embedding_stream(
        adata,
        color='condition',
        density=2,
        smooth=0.5,
        title='scVelo velocity',
        show=False,
    )
    plt.savefig("./outputs_of_scVelo/streamline_prediction_of_UniTVelo2.svg",
                bbox_inches='tight',
                pad_inches=0.1,
                transparent=True,
                dpi=300)
    plt.close()

    # 散点图
    scv.pl.scatter(
        adata,
        color='latent_time',
        color_map='plasma',
        size=20,
        title='scVelo latent time',
        show=False,
    )
    plt.savefig("./outputs_of_scVelo/umap_of_latent_time.svg",
                format='svg',
                bbox_inches='tight',
                pad_inches=0.1,
                transparent=True,
                dpi=300)
    plt.close()