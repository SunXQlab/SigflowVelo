import scanpy as sc
import scvelo as scv
import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def main():
    adata = sc.read("D:/TFvelo/TFvelo/data/TFvelo.h5ad")
    print(adata)

    # 你指定：velocity 在 layers['velocity']
    vkey = "velocity"

    if vkey not in adata.layers:
        raise KeyError(f"adata.layers 中没有 '{vkey}'，请确认 velocity layer 的名字。")

    # 避免 Windows 多进程问题
    scv.settings.n_jobs = 1

    # -----------------------------
    # 1) 计算 latent time（关键步骤）
    # -----------------------------
    # 这一步会在 adata.obs 中写入 'latent_time'
    scv.tl.latent_time(adata, vkey=vkey)

    lt_key = "latent_time"
    if lt_key not in adata.obs:
        raise RuntimeError("未在 adata.obs 中生成 'latent_time")

    print("latent_time summary:")
    print(adata.obs[lt_key].describe())

    # -----------------------------
    # 2) 确保 X_umap 可被 scVelo 正常切片
    # -----------------------------
    if "X_umap" not in adata.obsm:
        raise KeyError("adata.obsm 中没有 'X_umap'，无法 basis='umap' 作图。")

    if isinstance(adata.obsm["X_umap"], pd.DataFrame):
        adata.obsm["X_umap"] = adata.obsm["X_umap"].to_numpy()
    else:
        adata.obsm["X_umap"] = np.asarray(adata.obsm["X_umap"])

    # -----------------------------
    # 3) 可视化 latent time 并保存 SVG
    # -----------------------------
    ax = scv.pl.scatter(
        adata,
        basis="umap",
        color=lt_key,
        color_map="viridis",
        size=80,
        show=False
    )
    ax.set_title("Latent time of TFvelo")

    fig = ax.figure
    fig.canvas.draw()
    fig.savefig(
        "D:/TFvelo/TFvelo/results/TFvelo_latent_time.svg",
        format="svg",
        transparent=True
    )
    plt.close(fig)

    # -----------------------------
    # 4) 保存带 latent_time 的对象（可选）
    # -----------------------------
    adata.write("D:/TFvelo/TFvelo/data/TFvelo_with_latent_time.h5ad")

if __name__ == "__main__":
    mp.freeze_support()
    main()
