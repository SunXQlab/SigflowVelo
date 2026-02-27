import pandas as pd
import numpy as np
import os
import scanpy as sc
import scvelo as scv
import torch
import scipy
from scipy.spatial import distance_matrix
import Function

from sigflowvelo import SingleVelocity, calculate_velocity, plot_loss

adata_raw = sc.read('../output/dutcal/scPDAC_dutcal_final.h5ad')
gene_expression_matrix = adata_raw.X.toarray() if isinstance(adata_raw.X, (scipy.sparse.csr_matrix, scipy.sparse.csc_matrix)) else adata_raw.X
gene_names = adata_raw.var.index 
cell_names = adata_raw.obs.index  

df_count = pd.read_csv('../CCCInputData/Dutcal/count_file.csv', index_col=0)
df_imput = pd.read_csv('../CCCInputData/Dutcal/scImputeOutput/scimpute_count.csv', index_col=0)
InGem_link = pd.read_csv('../CCCInputData/Dutcal/Inflow-Gem_link_flie.csv', index_col=0)
GemOut_link = pd.read_csv('../CCCInputData/Dutcal/Gem-Outflow_link_flie.csv', index_col=0)

adata = sc.AnnData(X=df_count.values.astype(np.float64))  
adata.obs_names = df_count.index  
adata.var_names = df_count.columns  
adata.obs['celltype'] = adata_raw.obs["celltype"].copy()
adata.obs['group'] = adata_raw.obs["group"].copy()
adata.layers['Imputate'] = df_imput.values

Inflows = list(np.unique(InGem_link['Inflow'].values))
GEMs = list(np.unique(GemOut_link['GEM'].values))
Outflows = list(np.unique(GemOut_link['Outflow'].values))

n_gene = adata.shape[1]
adata.var['Inflows'] = np.full(n_gene, False, dtype=bool).astype(int)
adata.var['GEMs'] = np.full(n_gene, False, dtype=bool).astype(int)
adata.var['Outflows'] = np.full(n_gene, False, dtype=bool).astype(int)

for gene in list(adata.var_names):
    if gene in Inflows: adata.var['Inflows'][gene] = 1
    if gene in GEMs: adata.var['GEMs'][gene] = 1
    if gene in Outflows: adata.var['Outflows'][gene] = 1

adata.varm['GemOut_pair'] = np.full([n_gene, len(GEMs)], 'blank')
adata.varm['GemOut_regulate'] = np.full([n_gene, len(GEMs)], 0)

gene_names = list(adata.var_names)
for target in Outflows:
    if target in gene_names:
        target_idx = gene_names.index(target)
        df_tf_idx = np.where(GemOut_link['Outflow'].values == target)[0]
        tf_name = list(GemOut_link['GEM'].values[df_tf_idx])
        tf_idx = [index for index, element in enumerate(GEMs) if element in tf_name]

        for item1, item2 in zip(tf_idx, tf_name):
            adata.varm['GemOut_pair'][target_idx][item1] = item2
            adata.varm['GemOut_regulate'][target_idx][item1] = 1

InflowGem_score_flie = "../CCCInputData/Dutcal/Inflow-Gem_score_flie.csv"
InflowGem_score_matrix = pd.read_csv(InflowGem_score_flie, header=None, skiprows=1, usecols=range(1, 32)).values
InflowGem_allscore = np.tile(InflowGem_score_matrix, (len(adata.obs_names), 1, 1)) 
adata.obsm['InflowGem_signaling_score'] = InflowGem_allscore

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write_h5ad("../CCCInputData/Dutcal/flowsig_dutcal_CCCInput.h5ad")  

adata = sc.read('../CCCInputData/Dutcal/flowsig_dutcal_CCCInput.h5ad')
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

GEMs_expr = torch.tensor(adata.layers['Imputate'][:, adata.var['GEMs'].astype(bool)])
Outflows_expr = torch.tensor(adata.layers['Imputate'][:, adata.var['Outflows'].astype(bool)])

GemOut_regulate = torch.tensor(adata.varm['GemOut_regulate'])
nonzero_idx = torch.nonzero(torch.sum(GemOut_regulate, dim=1)).squeeze()
GemOut_regulate = GemOut_regulate[nonzero_idx].float()  

InGem_allscore = torch.tensor(adata.obsm['InflowGem_signaling_score'])

adata = Function.root_cell(adata, 0)
iroot = torch.tensor(adata.uns['iroot'])
print('The root cell is:', adata.uns['iroot'])

N_Outflows = Outflows_expr.shape[1]
layers_config = [N_Outflows+1, 32, 16, N_Outflows]

lr = 0.001
Weight = [0.1, 0.1, 0.0001, 0.001]

print('Initializing SingleVelocity model...')
model_Sen2Res = SingleVelocity(
    adata=adata,
    Outflows_expr=Outflows_expr,
    GEMs_expr=GEMs_expr,
    InGem_allscore=InGem_allscore,
    GemOut_regulate=GemOut_regulate,
    iroot=iroot,
    layers=layers_config,
    lr=lr,
    Weight=Weight,
    device=device
)

loss_history = model_Sen2Res.train(nIter=300, isSen2Res=True)

os.makedirs("../output/dutcal/", exist_ok=True)
torch.save(model_Sen2Res, "../output/dutcal/model_dutcal_Sen2Res_trained_v3.pth")

os.makedirs('../output/dutcal/Plot/', exist_ok=True)
plot_loss(loss_history, '../output/dutcal/Plot/LossAdam_Sen2Res.png')

adata_copy = adata.copy()
adata_copy.obsm['X_umap'] = adata_raw.obsm['X_umap']

print("Calculating velocity...")
adata_Sen2Res = calculate_velocity(adata_copy, model_Sen2Res, isSen2Res=True)
