import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42

import os
import glob
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import spatialdm as sdm
import spatialdm.plottings as pl

from spatialdm.plottings import plt_util

adata_ref = sc.read_h5ad("./ref_data/sc_cell2location.h5ad")
mod = cell2location.models.RegressionModel.load("./ref_data/", adata_ref)

adata_ref = mod.export_posterior(
    adata_ref, use_quantiles=True,
    add_to_varm=["q05","q50", "q95", "q0001"],
    sample_kwargs={'batch_size': 2500}
)

mod.plot_history(20)

if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:15]
inf_aver.to_csv("./Result/inf_aver.csv", index=True)

adata_ref.uns['mod']['factor_names']

def plot_genes_per_cell_type_embedding(slide, genes, ctypes):
    import scanpy as sc
    import matplotlib.pyplot as plt
    import numpy as np

    slide.var['SYMBOL'] = slide.var.index
    n_genes = len(genes)
    n_ctypes = len(ctypes)
    fig, axs = plt.subplots(
        nrows=n_genes, ncols=n_ctypes + 1,
        figsize=(4.5 * (n_ctypes + 1) + 2, 5 * n_genes + 1),
        squeeze=False
    )

    if 'spatial' not in slide.obsm:
        raise ValueError("slide.obsm['spatial'] not found. Please provide spatial coordinates.")

    slide.obsm['X_emb'] = slide.obsm['spatial']

    for j, gene in enumerate(genes):
        quantile_across_ct = np.array([
            np.quantile(slide.layers[n][:, slide.var["SYMBOL"] == gene].toarray(), 0.992)
            for n in slide.uns["mod"]["factor_names"]
        ])
        quantile_across_ct = np.partition(quantile_across_ct.flatten(), -2)[-2]

        sc.pl.embedding(
            slide,
            basis='X_emb',
            color=gene,
            gene_symbols='SYMBOL',
            size=20,
            cmap='magma',
            vmin=0,
            vmax='p99.2',
            ax=axs[j, 0],
            show=False
        )
        
        for i, ct in enumerate(ctypes):
            sc.pl.embedding(
                slide,
                basis='X_emb',
                color=gene,
                gene_symbols='SYMBOL',
                layer=ct,
                size=20,
                cmap='magma',
                vmin=0,
                vmax=quantile_across_ct,
                ax=axs[j, i + 1],
                show=False
            )
            axs[j, i + 1].set_title(f"{gene} {ct}")

    return fig, axs

adata_st = sc.read_h5ad("./cell2location_map/GSM8552941_PA02/st_cell2location_res.h5ad")
mod = cell2location.models.Cell2location.load("./cell2location_map/GSM8552941_PA02/", adata_st)

adata_st = mod.export_posterior(
    adata_st, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
)
mod.plot_QC()

adata_st.obsm

import pandas as pd
pd.DataFrame(adata_st.obsm['q05_cell_abundance_w_sf']).to_csv("./cell2location_map/st_cell2location_res.csv")

adata_st.obs[adata_st.uns['mod']['factor_names']] = adata_st.obsm['q05_cell_abundance_w_sf']
adata_st.obsm['q05_cell_abundance_w_sf']
adata_st.obs = adata_st.obs.join(adata_st.obsm["q05_cell_abundance_w_sf"])

from cell2location.utils import select_slide

adata_st
adata_st.obsm['X_emb'] = adata_st.obsm['spatial']

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.embedding(
        adata_st, basis='spatial',
        color=['Dutcal', 'Macrophage', 'T cells'],
        cmap='magma', size=20
    )
plt.savefig("./test_result/my_spatial.png", dpi=300, bbox_inches="tight")
plt.close()


mod.samples = adata_st.uns['mod']

expected_dict = mod.module.model.compute_expected_per_cell_type(
    mod.samples["post_sample_q05"], mod.adata_manager
)

for i, n in enumerate(mod.factor_names_):
    adata_st.layers[n] = expected_dict['mu'][i]

ctypes = ['Dutcal', 'Macrophage', 'T cells']
genes = ["SPP1", "CD44"]

with mpl.rc_context({'axes.facecolor':  'black'}):
    plot_genes_per_cell_type_embedding(adata_st, genes,ctypes);

plt.savefig("./test_result/gene_exp.png", dpi=300, bbox_inches="tight")
plt.close()

def plot_selected_pair(sample, pair, spots, selected_ind, figsize, cmap, cmap_l, cmap_r, 
                       celltype1=None, celltype2=None, title_fontsize=12, **kwargs):

    plot_kwargs = kwargs.copy()
    plot_kwargs.pop('title_fontsize', None)
    plot_kwargs.pop('celltype1', None)
    plot_kwargs.pop('celltype2', None)
    
    i = pd.Series(selected_ind == pair).idxmax()
    L = sample.uns['ligand'].loc[pair].dropna().values
    R = sample.uns['receptor'].loc[pair].dropna().values
    l1, l2 = len(L), len(R)
    
    if isinstance(sample.obsm['spatial'], pd.DataFrame):
        spatial_loc = sample.obsm['spatial'].values
    else:
        spatial_loc = sample.obsm['spatial']
    
    def extract_gene_from_pair(pair_name, gene_list, is_ligand=True):
        """
        Extract base gene name from a ligand–receptor pair label.
        """
        if '_' in pair_name:
            parts = pair_name.split('_')
            if is_ligand and len(parts) > 0:
                return parts[0]
            elif not is_ligand and len(parts) > 1:
                return parts[1]
            elif not is_ligand and len(parts) == 1:
                return parts[0]
        
        if len(gene_list) > 0:
            gene = gene_list[0]
            if ' ' in gene:
                return gene.split()[0]
            return gene
        
        return pair_name
    
    ligand_gene_name = extract_gene_from_pair(pair, L, is_ligand=True) if len(L) > 0 else "Unknown"
    receptor_gene_name = extract_gene_from_pair(pair, R, is_ligand=False) if len(R) > 0 else "Unknown"
    
    n_subplots = 1 + l1 + l2
    n_cols = min(n_subplots, 5)
    
    if isinstance(figsize, (tuple, list)) and len(figsize) == 2:
        subplot_width = figsize[0] / 5.0
        adjusted_width = subplot_width * n_cols
        adjusted_figsize = (adjusted_width, figsize[1])
    else:
        adjusted_figsize = figsize
    
    plt.figure(figsize=adjusted_figsize)
    plt.subplot(1, n_cols, 1)
    plt.scatter(spatial_loc[:,0], spatial_loc[:,1], c=spots.loc[pair], cmap=cmap,
                vmax=1, **plot_kwargs)
    plt_util('Moran: ' + str(sample.uns['local_stat']['n_spots'].loc[pair]) + ' spots')
    ax = plt.gca()
    ax.set_title(ax.get_title(), fontsize=title_fontsize)
    
    for l in range(l1):
        plt.subplot(1, n_cols, 2 + l)
        plt.scatter(spatial_loc[:,0], spatial_loc[:,1], c=sample[:,L[l]].X.toarray(),
                    cmap=cmap_l, **plot_kwargs)
        if celltype1:
            title_text = f'Ligand: {ligand_gene_name} in {celltype1}'
        else:
            title_text = 'Ligand: ' + L[l]
        plt_util(title_text)
        ax = plt.gca()
        ax.set_title(ax.get_title(), fontsize=title_fontsize)
    
    for l in range(l2):
        plt.subplot(1, n_cols, 2 + l1 + l)
        plt.scatter(spatial_loc[:,0], spatial_loc[:,1], c=sample[:,R[l]].X.toarray(),
                    cmap=cmap_r, **plot_kwargs)
        if celltype2:
            title_text = f'Receptor: {receptor_gene_name} in {celltype2}'
        else:
            title_text = 'Receptor: ' + R[l]
        plt_util(title_text)
        ax = plt.gca()
        ax.set_title(ax.get_title(), fontsize=title_fontsize)
    
def plot_pairs_fixed_2(sample, pairs_to_plot, pdf=None, figsize=(35, 5),
                     cmap='Greens', cmap_l='Reds', cmap_r='Reds', title='',
                     celltype1=None, celltype2=None, title_fontsize=12, **kwargs):

    if sample.uns['local_stat']['local_method'] == 'z-score':
        selected_ind = sample.uns['local_z_p'].index
        spots = 1 - sample.uns['local_z_p']
    if sample.uns['local_stat']['local_method'] == 'permutation':
        selected_ind = sample.uns['local_perm_p'].index
        spots = 1 - sample.uns['local_perm_p']

    if pdf is not None:
        with PdfPages(pdf + '.pdf') as pdf_writer:
            for pair in pairs_to_plot:
                plot_selected_pair(sample, pair, spots, selected_ind,
                                   figsize, cmap=cmap, cmap_l=cmap_l,
                                   cmap_r=cmap_r, celltype1=celltype1, 
                                   celltype2=celltype2, title_fontsize=title_fontsize, **kwargs)
                plt.tight_layout()
                pdf_writer.savefig(bbox_inches='tight')
                plt.close()
    else:
        for pair in pairs_to_plot:
            plot_selected_pair(sample, pair, spots, selected_ind,
                               figsize, cmap=cmap, cmap_l=cmap_l,
                               cmap_r=cmap_r, celltype1=celltype1,
                               celltype2=celltype2, title_fontsize=title_fontsize, **kwargs)
            plt.tight_layout()
            plt.show()
            plt.close()


def extract_lr_fixed(adata, species, mean='algebra', min_cell=0):
    """
    find overlapping LRs from local CellChatDB file
    """
    import pandas as pd
    from itertools import zip_longest
    import numpy as np

    if mean == 'geometric':
        from scipy.stats.mstats import gmean

    adata.uns['mean'] = mean

    lr_file = "/home/lsl/Bio/PDAC/cell2location/interaction_input_CellChatDB.csv"
    comp_file = "/home/lsl/Bio/PDAC/cell2location/complex_input_CellChatDB.csv"

    geneInter = pd.read_csv(lr_file, header=0, index_col=0)
    comp = pd.read_csv(comp_file, header=0, index_col=0)

    geneInter = geneInter.sort_values('annotation')
    ligand = geneInter.ligand.values
    receptor = geneInter.receptor.values
    geneInter.pop('ligand')
    geneInter.pop('receptor')

    t = []
    for i in range(len(ligand)):
        for n in [ligand, receptor]:
            l = n[i]
            if l in comp.index:
                n[i] = comp.loc[l].dropna().values[pd.Series(
                    comp.loc[l].dropna().values
                ).isin(adata.var_names)]
            else:
                n[i] = pd.Series(l).values[pd.Series(l).isin(adata.var_names)]

        if (len(ligand[i]) > 0) and (len(receptor[i]) > 0):
            if mean == 'geometric':
                meanL = gmean(adata[:, ligand[i]].X, axis=1)
                meanR = gmean(adata[:, receptor[i]].X, axis=1)
            else:
                meanL = adata[:, ligand[i]].X.mean(axis=1)
                meanR = adata[:, receptor[i]].X.mean(axis=1)

            if (sum(meanL > 0) >= min_cell) and (sum(meanR > 0) >= min_cell):
                t.append(True)
            else:
                t.append(False)
        else:
            t.append(False)

    ind = geneInter[t].index
    adata.uns['ligand'] = pd.DataFrame.from_records(
        zip_longest(*pd.Series(ligand[t]).values)
    ).transpose()
    adata.uns['ligand'].columns = [
        f'Ligand{i}' for i in range(adata.uns['ligand'].shape[1])
    ]
    adata.uns['ligand'].index = ind

    adata.uns['receptor'] = pd.DataFrame.from_records(
        zip_longest(*pd.Series(receptor[t]).values)
    ).transpose()
    adata.uns['receptor'].columns = [
        f'Receptor{i}' for i in range(adata.uns['receptor'].shape[1])
    ]
    adata.uns['receptor'].index = ind

    adata.uns['num_pairs'] = len(ind)
    adata.uns['geneInter'] = geneInter.loc[ind]

    if adata.uns['num_pairs'] == 0:
        raise ValueError("No effective RL pairs found. Please check input count matrix or file contents.")

    print(f"✅ 成功提取 {adata.uns['num_pairs']} 个有效配体-受体对")
    return

adata_M_D = adata_st.copy()
spp1_idx = np.where(adata_M_D.var['SYMBOL'] == 'SPP1')[0][0]
cd44_idx = np.where(adata_M_D.var['SYMBOL'] == 'CD44')[0][0]
adata_M_D.X[:, spp1_idx] = adata_M_D.layers['Macrophage'][:, spp1_idx].toarray().ravel()  
adata_M_D.X[:, cd44_idx] = adata_M_D.layers['Dutcal'][:, cd44_idx].toarray().ravel() 

adata_M_D.raw = adata_M_D


print("计算权重矩阵...")
sdm.weight_matrix(adata_M_D, l=1.2, cutoff=0.2, single_cell=False)

print("提取LR...")
sdm.extract_lr(adata_M_D, 'human', min_cell=3)

print("全局Moran selection...")
sdm.spatialdm_global(adata_M_D, 1000, specified_ind=None, method='both', nproc=1)

print("选择显著配对...")
sdm.sig_pairs(adata_M_D, method='permutation', fdr=True, threshold=0.1)

print("局部Spot selection...")
sdm.spatialdm_local(adata_M_D, n_perm=1000, method='both', specified_ind=None, nproc=1)

print("选择显著局部Spot...")
sdm.sig_spots(adata_M_D, method='permutation', fdr=False, threshold=0.1)

plot_pairs_fixed_2(adata_M_D, ['SPP1_CD44'], marker='s', pdf="local_pairs_PA02_SPP1_in_Macrophage_CD44_in_Ductal", 
                  celltype1="Macrophage", celltype2="Ductal")

adata_M_T = adata_st.copy()
spp1_idx = np.where(adata_M_T.var['SYMBOL'] == 'SPP1')[0][0]
cd44_idx = np.where(adata_M_T.var['SYMBOL'] == 'CD44')[0][0]
adata_M_T.X[:, spp1_idx] = adata_M_T.layers['Macrophage'][:, spp1_idx].toarray().ravel()  
adata_M_T.X[:, cd44_idx] = adata_M_T.layers['T cells'][:, cd44_idx].toarray().ravel() 

adata_M_T.raw = adata_M_T

print("计算权重矩阵...")
sdm.weight_matrix(adata_M_T, l=1.2, cutoff=0.2, single_cell=False)

print("提取LR...")
sdm.extract_lr(adata_M_T, 'human', min_cell=3)

print("全局Moran selection...")
sdm.spatialdm_global(adata_M_T, 1000, specified_ind=None, method='both', nproc=1)

print("选择显著配对...")
sdm.sig_pairs(adata_M_T, method='permutation', fdr=True, threshold=0.1)

print("局部Spot selection...")
sdm.spatialdm_local(adata_M_T, n_perm=1000, method='both', specified_ind=None, nproc=1)

print("选择显著局部Spot...")
sdm.sig_spots(adata_M_T, method='permutation', fdr=False, threshold=0.1)

plot_pairs_fixed_2(adata_M_T, ['SPP1_CD44'], marker='s', pdf="local_pairs_PA02_SPP1_in_Macrophage_CD44_in_T",
                  celltype1="Macrophage", celltype2="T cells")


adata_T_D = adata_st.copy()
spp1_idx = np.where(adata_T_D.var['SYMBOL'] == 'SPP1')[0][0]
cd44_idx = np.where(adata_T_D.var['SYMBOL'] == 'CD44')[0][0]
adata_T_D.X[:, spp1_idx] = adata_T_D.layers['T cells'][:, spp1_idx].toarray().ravel()  
adata_T_D.X[:, cd44_idx] = adata_T_D.layers['Dutcal'][:, cd44_idx].toarray().ravel() 

adata_T_D.raw = adata_T_D

print("计算权重矩阵...")
sdm.weight_matrix(adata_T_D, l=1.2, cutoff=0.2, single_cell=False)

print("提取LR...")
extract_lr_fixed(adata_T_D, 'human', min_cell=3)

print("全局Moran selection...")
sdm.spatialdm_global(adata_T_D, 1000, specified_ind=None, method='both', nproc=1)

print("选择显著配对...")
sdm.sig_pairs(adata_T_D, method='permutation', fdr=True, threshold=0.1)

print("局部Spot selection...")
sdm.spatialdm_local(adata_T_D, n_perm=1000, method='both', specified_ind=None, nproc=1)

print("选择显著局部Spot...")
sdm.sig_spots(adata_T_D, method='permutation', fdr=False, threshold=0.1)

plot_pairs_fixed_2(adata_T_D, ['SPP1_CD44'], marker='s', pdf="local_pairs_PA02_SPP1_in_T_CD44_in_Ductal",
                celltype1="T cells", celltype2="Ductal")


adata_T_M = adata_st.copy()
spp1_idx = np.where(adata_T_M.var['SYMBOL'] == 'SPP1')[0][0]
cd44_idx = np.where(adata_T_M.var['SYMBOL'] == 'CD44')[0][0]
adata_T_M.X[:, spp1_idx] = adata_T_M.layers['T cells'][:, spp1_idx].toarray().ravel()  
adata_T_M.X[:, cd44_idx] = adata_T_M.layers['Macrophage'][:, cd44_idx].toarray().ravel() 

adata_T_M.raw = adata_T_M

print("计算权重矩阵...")
sdm.weight_matrix(adata_T_M, l=1.2, cutoff=0.2, single_cell=False)

print("提取LR...")
extract_lr_fixed(adata_T_M, 'human', min_cell=3)

print("全局Moran selection...")
sdm.spatialdm_global(adata_T_M, 1000, specified_ind=None, method='both', nproc=1)

print("选择显著配对...")
sdm.sig_pairs(adata_T_M, method='permutation', fdr=True, threshold=0.1)

print("局部Spot selection...")
sdm.spatialdm_local(adata_T_M, n_perm=1000, method='both', specified_ind=None, nproc=1)

print("选择显著局部Spot...")
sdm.sig_spots(adata_T_M, method='permutation', fdr=False, threshold=0.1)

plot_pairs_fixed_2(adata_T_M, ['SPP1_CD44'], marker='s', pdf="local_pairs_PA02_SPP1_in_T_CD44_in_Macrophage",
                celltype1="T cells", celltype2="Macrophage")