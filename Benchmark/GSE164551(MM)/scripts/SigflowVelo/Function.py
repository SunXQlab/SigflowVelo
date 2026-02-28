def root_cell(adata, select_root):

    from scipy.spatial import distance_matrix
    import numpy as np
    import scanpy as sc
    import os
   
    # 简化：只从"condition"为"sensitive"的细胞中选择root cell
    if 'condition' in adata.obs.columns:
        # 创建只包含sensitive细胞的子集
        sensitive_mask = adata.obs['condition'] == 'sensitive'
        adata_sensitive = adata[sensitive_mask].copy()
        
        if adata_sensitive.shape[0] == 0:
            raise ValueError("No cells with condition 'sensitive' found in the dataset")
        
        print(f"Using {adata_sensitive.shape[0]} cells with condition 'sensitive' for root cell selection")
    else:
        raise ValueError("'condition' column not found in adata.obs")
   
    # 后续所有操作都在adata_sensitive子集上进行
    if select_root == 'STAGATE':
        max_cell_for_subsampling = 5000
        if adata_sensitive.shape[0] < max_cell_for_subsampling:
            sub_adata_x = adata_sensitive.obsm['X_umap']
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            sensitive_root_index = np.argmax(sum_dists)
        else:
            indices = np.arange(adata_sensitive.shape[0])
            selected_ind = np.random.choice(indices, max_cell_for_subsampling, False)
            sub_adata_x = adata_sensitive.obsm['X_umap'][selected_ind, :]
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            sensitive_root_index = np.argmax(sum_dists)

    elif select_root == 'UMAP':
        adata_copy = adata_sensitive.copy()

        keys_to_delete = ['X_umap', 'neighbors', 'umap', 'distances', 'connectivities']
        for key in keys_to_delete:
            if key in adata_copy.obsm:
                del adata_copy.obsm[key]
            if key in adata_copy.uns:
                del adata_copy.uns[key]
            if key in adata_copy.obsp:
                del adata_copy.obsp[key]

        sc.tl.pca(adata_copy, svd_solver="arpack")
        sc.pp.neighbors(adata_copy, n_pcs=50)
        sc.tl.umap(adata_copy)

        max_cell_for_subsampling = 5000
        if adata_copy.shape[0] < max_cell_for_subsampling:
            sub_adata_x = adata_copy.obsm['X_umap']
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            sensitive_root_index = np.argmax(sum_dists)
        else:
            indices = np.arange(adata_copy.shape[0])
            selected_ind = np.random.choice(indices, max_cell_for_subsampling, False)
            sub_adata_x = adata_copy.obsm['X_umap'][selected_ind, :]
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            sensitive_root_index = np.argmax(sum_dists)
            
    elif select_root == 'CCC_genes':
        lig = adata_sensitive.var['ligand'].astype(bool)
        rec = adata_sensitive.var['receptor'].astype(bool)
        tf = adata_sensitive.var['TFs'].astype(bool)
        tg = adata_sensitive.var['TGs'].astype(bool)
        combined_bool = lig | rec | tf | tg
        sub_adata = adata_sensitive[:, combined_bool]
        sub_adata = np.unique(sub_adata.var_names)
        sub_adata = adata_sensitive[:, sub_adata]
        max_cell_for_subsampling = 5000
        if sub_adata.shape[0] < max_cell_for_subsampling:
            sub_adata_x = sub_adata.layers['Imputate']
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            sensitive_root_index = np.argmax(sum_dists)
        else:
            indices = np.arange(sub_adata.shape[0])
            selected_ind = np.random.choice(indices, max_cell_for_subsampling, False)
            sub_adata_x = sub_adata.layers['Imputate'][selected_ind, :]
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            sensitive_root_index = np.argmax(sum_dists)
            
    elif select_root == "spatial":
        max_cell_for_subsampling = 50000
        if adata_sensitive.shape[0] < max_cell_for_subsampling:
            sub_adata_x = adata_sensitive.obsm['spatial']
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            sensitive_root_index = np.argmax(sum_dists)
        else:
            indices = np.arange(adata_sensitive.shape[0])
            selected_ind = np.random.choice(indices, max_cell_for_subsampling, False)
            sub_adata_x = adata_sensitive.obsm['spatial'][selected_ind, :]
            sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
            sensitive_root_index = np.argmax(sum_dists)
            
    elif type(select_root) == int:
        # 检查指定的root cell是否在sensitive细胞中
        if select_root >= adata_sensitive.shape[0]:
            raise ValueError(f"The specified root cell index {select_root} is out of range for sensitive cells")
        sensitive_root_index = select_root

    # 将sensitive子集中的索引转换为原始adata中的索引
    sensitive_indices = np.where(sensitive_mask)[0]
    adata.uns['iroot'] = sensitive_indices[sensitive_root_index]
    print(f"Selected root cell index: {adata.uns['iroot']}")
    
    return adata

def create_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Directory created at: {path}")

def save_model_and_data(model, data, path):
    torch.save(model, os.path.join(path, "model_spa_velo.pth"))
    torch.save(data, os.path.join(path, "CCCvelo.pt"))
    # print(f"Model and data saved at: {path}")