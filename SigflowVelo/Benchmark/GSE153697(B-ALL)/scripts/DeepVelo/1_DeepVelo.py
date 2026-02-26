import deepvelo as dv
import scvelo as scv
from deepvelo.utils.temporal import latent_time

if __name__ == "__main__":
    adata = scv.read("../../data/scPDAC_with_su_normalized.h5ad")
    print (adata)
    scv.pp.moments(adata, n_neighbors=30, n_pcs=30)
    # run DeepVelo using the default configs
    trainer = dv.train(adata, dv.Constants.default_configs)
    adata=(latent_time(adata,copy=True))
    # this will train the model and predict the velocity vectore. The result is stored in adata.layers['velocity']. You can use trainer.model to access the model.
    adata.write_h5ad("../../results/DeepVelo/adata/adata_with_DeepVelo.h5ad")
    print (adata)