import scanpy as sc
import numpy as np
from scipy.spatial import distance_matrix 
import os
import torch

def root_cell(adata, select_root):

    if select_root == 'STAGATE':
        max_cell_for_subsampling = 5000
        if adata.shape[0] < max_cell_for_subsampling:
            sub_adata_x = adata.obsm['X_umap']
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            adata.uns['iroot'] = np.argmax(sum_dists)
        else:
            indices = np.arange(adata.shape[0])
            selected_ind = np.random.choice(indices, max_cell_for_subsampling, False)
            sub_adata_x = adata.obsm['X_umap'][selected_ind, :]
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            adata.uns['iroot'] = np.argmax(sum_dists)
    elif select_root == 'UMAP':
        adata_copy = adata.copy()

# --- 建议修改：安全地删除键 ---
        # 逐个尝试删除，如果键不存在 (KeyError) 则跳过
        try:
            del adata_copy.obsm['X_umap']
        except KeyError:
            pass  # 如果 'X_umap' 不存在，则忽略
        
        try:
            del adata_copy.uns['neighbors']
        except KeyError:
            pass  # 如果 'neighbors' 不存在，则忽略

        try:
            del adata_copy.uns['umap']
        except KeyError:
            pass  # 如果 'umap' 不存在，则忽略

        try:
            del adata_copy.obsp['distances']
        except KeyError:
            pass  # 如果 'distances' 不存在，则忽略

        try:
            del adata_copy.obsp['connectivities']
        except KeyError:
            pass  # 如果 'connectivities' 不存在，则忽略

        sc.tl.pca(adata_copy, svd_solver="arpack")
        sc.pp.neighbors(adata_copy, n_pcs=50)
        sc.tl.umap(adata_copy)

        max_cell_for_subsampling = 5000
        if adata_copy.shape[0] < max_cell_for_subsampling:
            sub_adata_x = adata_copy.obsm['X_umap']
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            adata.uns['iroot'] = np.argmax(sum_dists)
        else:
            indices = np.arange(adata_copy.shape[0])
            selected_ind = np.random.choice(indices, max_cell_for_subsampling, False)
            sub_adata_x = adata_copy.obsm['X_umap'][selected_ind, :]
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            adata.uns['iroot'] = np.argmax(sum_dists)
    elif select_root == 'CCC_genes':
        lig = adata.var['ligand'].astype(bool)
        rec = adata.var['receptor'].astype(bool)
        tf = adata.var['TFs'].astype(bool)
        tg = adata.var['TGs'].astype(bool)
        combined_bool = lig | rec | tf | tg
        sub_adata = adata[:, combined_bool]
        sub_adata = np.unique(sub_adata.var_names)
        sub_adata = adata[:, sub_adata]
        max_cell_for_subsampling = 5000
        if sub_adata.shape[0] < max_cell_for_subsampling:
            sub_adata_x = sub_adata.layers['Imputate']
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            adata.uns['iroot'] = np.argmax(sum_dists)
        else:
            indices = np.arange(sub_adata.shape[0])
            selected_ind = np.random.choice(indices, max_cell_for_subsampling, False)
            sub_adata_x = sub_adata.layers['Imputate'][selected_ind, :]
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            adata.uns['iroot'] = np.argmax(sum_dists)
    elif select_root == "spatial":
        max_cell_for_subsampling = 50000
        if adata.shape[0] < max_cell_for_subsampling:
            sub_adata_x = adata.obsm['spatial']
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            adata.uns['iroot'] = np.argmax(sum_dists)
        else:
            indices = np.arange(adata.shape[0])
            selected_ind = np.random.choice(indices, max_cell_for_subsampling, False)
            sub_adata_x = adata.obsm['spatial'][selected_ind, :]
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            adata.uns['iroot'] = np.argmax(sum_dists)
    elif type(select_root) == int:
        adata.uns['iroot'] = select_root

    return adata

def create_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Directory created at: {path}")

def save_model_and_data(model, data, path):
    torch.save(model, os.path.join(path, "model_spa_velo.pth"))
    torch.save(data, os.path.join(path, "CCCvelo.pt"))
    # print(f"Model and data saved at: {path}")
