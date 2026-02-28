import TFvelo as TFv
import anndata as ad
import scanpy as sc
import numpy as np
import scipy
import matplotlib
import pandas as pd
matplotlib.use('AGG')


np.set_printoptions(suppress=True)


def check_data_type(adata):
    import pandas as pd
    for key in list(adata.var):
        # Get unique values (convert to string for comparison)
        unique_values = set(str(val) for val in adata.var[key].unique())
        
        # Check if only contains 'True' and 'False'
        if unique_values.issubset({'True', 'False', 'nan'}):
            print('Checking', key)
            # Use safer conversion method
            def convert_value(val):
                if pd.isna(val):
                    return np.nan
                elif str(val) == 'True':
                    return True
                elif str(val) == 'False':
                    return False
                else:
                    return val
            
            adata.var[key] = adata.var[key].apply(convert_value)
    
    return adata


def data_type_tostr(adata, key=None):
    if key is None:
        for key in list(adata.var):
            if adata.var[key][0] in [True, False]:
                print('Transfering', key)
                adata.var[key] = adata.var[key].map({True: 'True', False:'False'})
    elif key in adata.var.keys():
        if adata.var[key][0] in [True, False]:
            print('Transfering', key)
            adata.var[key] = adata.var[key].map({True: 'True', False:'False'})
    return   



def get_pseudotime(adata):
    TFv.tl.velocity_graph(adata, basis=None, vkey='velocity', xkey='M_total')
    TFv.tl.velocity_pseudotime(adata, vkey='velocity', modality='M_total') 
    TFv.pl.scatter(adata, basis=args.basis, color='velocity_pseudotime', cmap='gnuplot', fontsize=20, save='pseudotime')
    return adata


def get_sort_positions(arr):
    positions = np.argsort(np.argsort(arr))
    positions_normed = positions/(len(arr)-1)
    return positions_normed


def get_metric_pseudotime(adata, t_key='latent_t'):
    n_cells, n_genes = adata.shape
    adata.var['spearmanr_pseudotime'] = 0.0
    for i in range(n_genes):
        correlation, _ = scipy.stats.spearmanr(adata.layers[t_key][:,i], adata.obs['velocity_pseudotime'])
        adata.var['spearmanr_pseudotime'][i] = correlation
    return adata

def show_adata(args, adata, save_name, show_all=0):
    """
    Only plot stream plot, not single gene plots
    This avoids categorical variable issues
    """
    
    # Ensure clusters categories are strings
    if 'clusters' in adata.obs and hasattr(adata.obs['clusters'], 'cat'):
        clusters = adata.obs['clusters']
        categories = clusters.cat.categories
        if len(categories) > 0 and not isinstance(categories[0], str):
            print("Converting clusters categories to strings")
            new_categories = categories.astype(str)
            adata.obs['clusters'] = clusters.cat.rename_categories(new_categories)
    
    # Ensure embedding coordinates are numpy arrays, not pandas DataFrame
    basis_key = f'X_{args.basis}'
    if basis_key in adata.obsm:
        import pandas as pd
        if isinstance(adata.obsm[basis_key], pd.DataFrame):
            print(f"Converting {basis_key} from DataFrame to numpy array")
            adata.obsm[basis_key] = adata.obsm[basis_key].values
        elif not isinstance(adata.obsm[basis_key], np.ndarray):
            print(f"Converting {basis_key} to numpy array")
            adata.obsm[basis_key] = np.array(adata.obsm[basis_key])
    
    # Only plot stream plot, skip single gene plots
    if len(adata.obs['clusters'].cat.categories) > 10:
        legend_loc = 'right margin'
    else:
        legend_loc = 'on data'
    
    cutoff_perc = 20
    
    try:
        print("Plotting stream plot...")
        TFv.pl.velocity_embedding_stream(
            adata, 
            vkey='velocity', 
            use_derivative=False, 
            density=4, 
            basis=args.basis,
            cutoff_perc=cutoff_perc, 
            smooth=0.5, 
            fontsize=5, 
            recompute=True,
            legend_loc=legend_loc, 
            color='condition',
            save='embedding_stream_' + save_name
        )
        print("Stream plot completed")
    except Exception as e:
        print(f"Error plotting stream plot: {e}")
        import traceback
        traceback.print_exc()
    
    return

def get_sort_t(adata):
    t = adata.layers['fit_t_raw'].copy()
    normed_t = adata.layers['fit_t_raw'].copy()
    n_bins = 20
    n_cells, n_genes = adata.shape
    sort_t = np.zeros([n_cells, n_genes])
    non_blank_gene = np.zeros(n_genes, dtype=int)
    hist_all, bins_all = np.zeros([n_genes, n_bins]), np.zeros([n_genes, n_bins+1])
    for i in range(n_genes):
        gene_name = adata.var_names[i]
        tmp = t[:,i].copy()
        if np.isnan(tmp).sum():
            non_blank_gene[i] = 1 
            continue
        hist, bins = np.histogram(tmp, bins=n_bins)
        hist_all[i], bins_all[i] = hist, bins
        if not (0 in list(hist)):
            if (tmp.min() < 0.1) and (tmp.max() > 0.8):
                blank_start_bin_id = np.argmin(hist)
                blank_end_bin_id = blank_start_bin_id
                non_blank_gene[i] = 1
                blank_start_bin = bins[blank_start_bin_id]
                blank_end_bin = bins[blank_end_bin_id]
                tmp = (tmp < blank_start_bin)*1 + tmp 
            else:
                blank_end_bin = tmp.min()
        else:
            blank_start_bin_id = list(hist).index(0)
            for j in range(blank_start_bin_id+1, len(hist)):
                if hist[j] > 0:
                    blank_end_bin_id = j
                    break
            blank_start_bin = bins[blank_start_bin_id]
            blank_end_bin = bins[blank_end_bin_id]
            tmp = (tmp < blank_start_bin)*1 + tmp 
            
        t[:,i] = tmp
        tmp = tmp - blank_end_bin
        tmp = tmp/tmp.max()
        normed_t[:,i] = tmp
        sort_t[:,i] = get_sort_positions(tmp)

    adata.layers['latent_t'] = sort_t.copy() 
    adata.var['non_blank_gene'] = non_blank_gene.copy()
    return adata


def main(args):
    # Read data from the path generated by the previous script
    adata = ad.read_h5ad(args.data_path+"rc.h5ad") 
    adata.var_names_make_unique()
    print("Output adata obtained from run_demo")
    print(adata)
    
    # Check if key data exists
    print("Checking if key data exists:")
    print(f"  'velo_hat' in layers: {'velo_hat' in adata.layers}")
    print(f"  'fit_scaling_y' in var: {'fit_scaling_y' in adata.var.columns}")
    
    if 'fit_scaling_y' in adata.var.columns:
        print(f"  fit_scaling_y dtype: {adata.var['fit_scaling_y'].dtype}")
        print(f"  fit_scaling_y first value: {adata.var['fit_scaling_y'].iloc[0]}")
    
    check_data_type(adata)
    print(adata)

    # Check key data again after conversion
    print("\nChecking key data after conversion:")
    print(f"  'velo_hat' in layers: {'velo_hat' in adata.layers}")
    print(f"  'fit_scaling_y' in var: {'fit_scaling_y' in adata.var.columns}")
    
    if 'fit_scaling_y' in adata.var.columns:
        print(f"  fit_scaling_y dtype: {adata.var['fit_scaling_y'].dtype}")
        print(f"  fit_scaling_y first value: {adata.var['fit_scaling_y'].iloc[0]}")

    # Create velocity layer
    import numpy as np
    if 'velocity' not in adata.layers:
        print("Creating 'velocity' layer...")
        # Ensure velo_hat and fit_scaling_y exist
        if 'velo_hat' in adata.layers and 'fit_scaling_y' in adata.var:
            # Ensure fit_scaling_y is numeric type
            try:
                fit_scaling_y_values = np.array(adata.var['fit_scaling_y'], dtype=float)
                n_cells = adata.shape[0]
                expanded_scaling_y = np.expand_dims(fit_scaling_y_values, 0).repeat(n_cells, axis=0)
                adata.layers['velocity'] = adata.layers['velo_hat'] / expanded_scaling_y
                print("'velocity' layer created successfully, shape:", adata.layers['velocity'].shape)
            except Exception as e:
                print(f"Error creating velocity layer: {e}")
                print(f"fit_scaling_y type: {type(adata.var['fit_scaling_y'])}")
                print(f"fit_scaling_y first 5 values: {adata.var['fit_scaling_y'].head()}")
                return
        else:
            print("Error: Missing 'velo_hat' or 'fit_scaling_y', cannot create velocity.")
            print(f"'velo_hat' in layers: {'velo_hat' in adata.layers}")
            print(f"'fit_scaling_y' in var: {'fit_scaling_y' in adata.var.columns}")
            return

    losses = adata.varm['loss'].copy()
    losses[np.isnan(losses)] = 1e6
    adata.var['min_loss'] = losses.min(1)

    if 'X_pca' not in adata.obsm.keys():
        print('Running PCA...')
        sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
    
    if (args.basis == 'umap') and ('X_umap' not in adata.obsm.keys()):
        print('Running UMAP...')
        if args.dataset_name == 'hesc1':
            sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
            sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=30, n_pcs=5)
            sc.tl.umap(adata)
        else:
            sc.tl.umap(adata)  
            sc.pl.umap(adata, color='clusters', save=True)
    
    adata = get_pseudotime(adata)
    
    adata_copy = adata.copy()
    adata_copy = get_sort_t(adata_copy) 

    TFv.tl.velocity_graph(adata_copy, basis=None, vkey='velo_hat', xkey='M_total')
    adata_copy.uns['clusters_colors'] = adata.uns['clusters_colors']
    show_adata(args, adata_copy, save_name='velo', show_all=0)

    # data_type_tostr(adata_copy)
    print(adata_copy)
    adata_copy.write(args.data_path + 'TFvelo.h5ad')

    return

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset_name', type=str, default='Multi_myeloma_NK', 
                       help='pancreas, gastrulation_erythroid, 10x_mouse_brain, hesc1') 
    parser.add_argument('--layer', type=str, default="M_total", help='M_total, total') 
    parser.add_argument('--basis', type=str, default="umap", help='umap, tsne, pca')
    parser.add_argument('--loss_percent_thres', type=int, default=50, help='max loss of each gene')
    parser.add_argument('--spearmanr_thres', type=float, default=0.8, help='min spearmanr')
    parser.add_argument('--save_name', type=str, default='_demo', help='save_name')
    args = parser.parse_args() 
    
    # Use the existing path rule for data saving
    args.data_path = 'TFvelo_'+ args.dataset_name + args.save_name+ '/'
    print('------------------------------------------------------------')   

    print(args) 
    main(args)