import scanpy as sc
import pandas as pd
import flowsig as fs

adata = sc.read('scPDAC_dutcal_normalized.h5ad')

condition_key = 'group'

cellchat_Sensitive = pd.read_csv('../CellChat/Output/scPDAC_celltype_communications_Sensitive.csv')
cellchat_Resist = pd.read_csv('../CellChat/Output/scPDAC_celltype_communications_Resist.csv')
cellchat_Sensitive = cellchat_Sensitive[(cellchat_Sensitive['source'] == 'Dutcal') | (cellchat_Sensitive['target'] == 'Dutcal')]
cellchat_Resist = cellchat_Resist[(cellchat_Resist['source'] == 'Dutcal') | (cellchat_Resist['target'] == 'Dutcal')]
cellchat_Sensitive.to_csv('../output/dutcal/cellchat_dutcal_sensitive.csv')
cellchat_Resist.to_csv('../output/dutcal/cellchat_dutcal_resist.csv')

cellchat_output_key = 'cellchat_output'

adata.uns[cellchat_output_key] = {'Sensitive': cellchat_Sensitive,
                                  'Resist': cellchat_Resist}

counts_matrix = adata.X.copy() 
adata.layers["counts"] = counts_matrix

if adata.obs.index.name is None:
    adata.obs.index.name = 'cell_id' 

if adata.var.index.name is None:
    adata.var.index.name = 'gene_id' 

fs.pp.construct_gems_using_pyliger(adata,
                                n_gems = 10,
                                layer_key = 'counts',
                                condition_key = condition_key)

gems_result = adata.obsm['X_gem']
gems_df = pd.DataFrame(gems_result, index=adata.obs_names, columns=[f"Module_{i+1}" for i in range(gems_result.shape[1])])
gems_df.to_csv("../output/dutcal/gems_result.csv")
adata.write_h5ad("../output/dutcal/scPDAC_dutcal+GEMs.h5ad")  

fs.pp.construct_flows_from_cellchat(adata,
                                cellchat_output_key,
                                gem_expr_key = 'X_gem',
                                scale_gem_expr = True,
                                model_organism = 'human',
                                flowsig_network_key = 'flowsig_network',
                                flowsig_expr_key = 'X_flow')

fs.pp.determine_informative_variables(adata,  
                                    flowsig_expr_key = 'X_flow',
                                    flowsig_network_key = 'flowsig_network',
                                    spatial = False,
                                    condition_key = condition_key,
                                    qval_threshold = 0.05,
                                    logfc_threshold = 0.5)

flow_var_info = adata.uns['flowsig_network']['flow_var_info']
flow_var_info_df = pd.DataFrame(flow_var_info)
flow_var_info_df.to_csv('../output/dutcal/flow_var_info.csv', index=True)

adata.write_h5ad("../output/dutcal/scPDAC_dutcal+GEMs+Signal.h5ad")  

fs.tl.learn_intercellular_flows(adata,
                        condition_key = condition_key,
                        control_key = 'Sensitive', 
                        flowsig_key = 'flowsig_network',
                        flow_expr_key = 'X_flow',
                        use_spatial = False,
                        n_jobs = 1,
                        n_bootstraps = 10)


fs.tl.apply_biological_flow(adata,
                        flowsig_network_key = 'flowsig_network',
                        adjacency_key = 'adjacency')

edge_threshold = 0.7

fs.tl.filter_low_confidence_edges(adata,
                                edge_threshold = edge_threshold,
                                flowsig_network_key = 'flowsig_network',
                                adjacency_key = 'adjacency')

adata.write('../output/dutcal/burkhardt21_merged.h5ad', compression='gzip')

adjacency = adata.uns['flowsig_network']['network']['adjacency']
adjacency_df = pd.DataFrame(adjacency, index=flow_var_info_df.index, columns=flow_var_info_df.index)
adjacency_df.to_csv('../output/dutcal/adjacency_network.csv', index=True)

adata.write_h5ad("../output/dutcal/scPDAC_dutcal_final.h5ad")  