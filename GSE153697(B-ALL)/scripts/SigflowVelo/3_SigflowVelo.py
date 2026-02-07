import scipy
import torch
from collections import OrderedDict
import scanpy as sc
import pandas as pd
# import TFvelo as TFv
import numpy as np
import random
import matplotlib.pyplot as plt
import warnings
import time
from tqdm import tqdm
from collections import Counter
import scvelo as scv
import seaborn as sns
import Function

# Set random seed for reproducibility
def set_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

class DNN(torch.nn.Module):
    def __init__(self, layers):
        super(DNN, self).__init__()

        # parameters
        self.depth = len(layers) - 1

        # set up layer order dict
        self.activation = torch.nn.Tanh

        layer_list = list()
        for i in range(self.depth - 1):
            layer_list.append(
                ('layer_%d' % i, torch.nn.Linear(layers[i], layers[i + 1]))
            )
            layer_list.append(('activation_%d' % i, self.activation()))

        layer_list.append(
            ('layer_%d' % (self.depth - 1), torch.nn.Linear(layers[-2], layers[-1]))
        )
        layerDict = OrderedDict(layer_list)

        # deploy layers
        self.layers = torch.nn.Sequential(layerDict)

    def forward(self, x):
        out = self.layers(x)
        return out

# Here TGs_expr represents Outflows_expr, TFs_expr represents GEMs_expr,
# TFLR_allscore represents InGem_allscore, TGTF_regulate represents GemOut_regulate
class SingleVelocity_v2():
    def __init__(self, TGs_expr, TFs_expr, TFLR_allscore, TGTF_regulate, iroot, layers, lr, Weight, device, adata):
        # data
        self.adata = adata
        self.TGs_expr = TGs_expr.clone().detach().float().to(device)
        self.TFs_expr = TFs_expr.clone().detach().float().to(device)
        self.TFLR_allscore = TFLR_allscore.clone().detach().float().to(device)
        self.regulate = TGTF_regulate.clone().detach().to(device)  # torch.Size([539, 114])
        self.iroot = iroot.int().to(device)
        # self.t = torch.linspace(0, 1, 2000).unsqueeze(1).requires_grad_(True).to(device)
        self.t1 = torch.clamp(torch.normal(mean = 0.75, std = 1, size=(1000,)), min=0.3, max = 1).unsqueeze(1).requires_grad_(True).to(device)
        self.t2 = torch.clamp(torch.normal(mean = 0.25, std = 1, size=(1000,)), min=0, max = 0.7).unsqueeze(1).requires_grad_(True).to(device)
        self.t = torch.cat([self.t1, self.t2], dim=0)
        self.Weight = Weight
        self.N_cell = TGs_expr.shape[0]
        self.N_TFs = TFs_expr.shape[1]
        self.N_TGs = TGs_expr.shape[1]
        self.N_LRs = TFLR_allscore.shape[2]

        self.rootcell_exp = self.TGs_expr[self.iroot, :]

        # setting parameters
        # self.V1 = torch.empty((self.N_TFs, self.N_LRs), dtype=torch.float32).uniform_(1, 10).float().requires_grad_(True).to(device)
        # self.K1 = torch.empty((self.N_TFs, self.N_LRs), dtype=torch.float32).uniform_(1, 10).float().requires_grad_(True).to(device)
        # self.V2 = torch.empty((self.N_TGs, self.N_TFs), dtype=torch.float32).uniform_(1, 10).float().requires_grad_(True).to(device)
        # self.K2 = torch.empty((self.N_TGs, self.N_TFs), dtype=torch.float32).uniform_(1, 10).float().requires_grad_(True).to(device)

        self.V1 = torch.empty((self.N_TFs, self.N_LRs), dtype=torch.float32).uniform_(0, 1).float().requires_grad_(True).to(device)
        self.K1 = torch.empty((self.N_TFs, self.N_LRs), dtype=torch.float32).uniform_(0, 1).float().requires_grad_(True).to(device)
        self.V2 = torch.empty((self.N_TGs, self.N_TFs), dtype=torch.float32).uniform_(0, 1).float().requires_grad_(True).to(device)
        self.K2 = torch.empty((self.N_TGs, self.N_TFs), dtype=torch.float32).uniform_(0, 1).float().requires_grad_(True).to(device)

        # self.V1 = torch.empty((self.N_TFs, self.N_LRs), dtype=torch.float32).uniform_(-1, 1).float().requires_grad_(True).to(device)
        # self.K1 = torch.empty((self.N_TFs, self.N_LRs), dtype=torch.float32).uniform_(-1, 1).float().requires_grad_(True).to(device)
        # self.V2 = torch.empty((self.N_TGs, self.N_TFs), dtype=torch.float32).uniform_(-1, 1).float().requires_grad_(True).to(device)
        # self.K2 = torch.empty((self.N_TGs, self.N_TFs), dtype=torch.float32).uniform_(-1, 1).float().requires_grad_(True).to(device)

        self.V1 = torch.nn.Parameter(self.V1)  # Backpropagation requires updated parameters
        self.K1 = torch.nn.Parameter(self.K1)
        self.V2 = torch.nn.Parameter(self.V2)
        self.K2 = torch.nn.Parameter(self.K2)

        # deep neural networks
        self.dnn = DNN(layers).to(device)
        self.dnn.register_parameter('V1', self.V1)
        self.dnn.register_parameter('K1', self.K1)
        self.dnn.register_parameter('V2', self.V2)
        self.dnn.register_parameter('K2', self.K2)

        self.optimizer_Adam = torch.optim.Adam(self.dnn.parameters(), lr=lr)
        self.iter = 0

    def net_z(self):
        t = self.t
        N_TGs = self.N_TGs
        z0 = self.rootcell_exp.repeat(t.size(0), 1)
        z_and_t = torch.cat([z0, t], dim=1)
        z_dnn = self.dnn(z_and_t)

        for i in range(N_TGs):
            z_t_pre = torch.autograd.grad(
                z_dnn[:, i], t,
                grad_outputs=torch.ones_like(z_dnn[:, i]),
                retain_graph=True,
                create_graph=True
            )[0]
            if i == 0:
                dz_dt = z_t_pre
            else:
                dz_dt = torch.cat((dz_dt, z_t_pre), 1)
        z_dnn = torch.where(z_dnn > 0, z_dnn, torch.full_like(z_dnn, 0))

        return z_dnn, dz_dt
    
    def net_z_v2(self, t):
        z0 = self.rootcell_exp.repeat(t.size(0), 1)
        z_and_t = torch.cat([z0, t], dim=1)
        z_dnn = self.dnn(z_and_t)
        z_dnn = torch.where(z_dnn > 0, z_dnn, torch.full_like(z_dnn, 0))

        return z_dnn
    
    def assign_latenttime_v3(self):
        tpoints_sens = self.t2  # t2（0.25）applied to "sensitive" cell
        tpoints_resist = self.t1  # t1（0.75） applied to "resistant" cell

        z_dnn_sens = self.net_z_v2(tpoints_sens)
        z_dnn_resist = self.net_z_v2(tpoints_resist)
        z_obs = self.TGs_expr

        # obtain the group information
        groups = self.adata.obs["group"]

        # assign latent time for sensitive cell
        sensitive_mask = (groups == "sensitive").values
        # print('the sensitive mask is:\n', sensitive_mask)
        if sensitive_mask.sum() > 0:
            z_obs_sens = z_obs[sensitive_mask]
            loss_cell_to_t_sens = torch.sum((z_dnn_sens.unsqueeze(1) - z_obs_sens.unsqueeze(0)) ** 2, dim=2)
            pos_sens = torch.argmin(loss_cell_to_t_sens, dim=0)
            pos_sens = torch.clamp(pos_sens, 0, tpoints_sens.size(0) - 1)
            fit_t_sens = tpoints_sens[pos_sens]
            fit_t_sens = fit_t_sens.flatten()[:, None].squeeze()
        else:
            fit_t_sens = None
            pos_sens = None

        # assign latent time for resistant cell
        resist_mask = (groups == "resistant").values
        # print('the resist mask is:\n', resist_mask)
        if resist_mask.sum() > 0:
            z_obs_resist = z_obs[resist_mask]
            loss_cell_to_t_resist = torch.sum((z_dnn_resist.unsqueeze(1) - z_obs_resist.unsqueeze(0)) ** 2, dim=2)
            pos_resist = torch.argmin(loss_cell_to_t_resist, dim=0)
            pos_resist = torch.clamp(pos_resist, 0, tpoints_resist.size(0) - 1)
            fit_t_resist = tpoints_resist[pos_resist]
            fit_t_resist = fit_t_resist.flatten()[:, None].squeeze()
        else:
            fit_t_resist = None
            pos_resist = None

        fit_t = torch.zeros(z_obs.shape[0], dtype=tpoints_sens.dtype)
        pos = torch.zeros(z_obs.shape[0], dtype=torch.long)

        if fit_t_sens is not None:
            fit_t[sensitive_mask] = fit_t_sens
            pos[sensitive_mask] = pos_sens
        if fit_t_resist is not None:
            fit_t[resist_mask] = fit_t_resist
            pos[resist_mask] = pos_resist
            
        # print('the sensitive fit latent time shape is:\n', fit_t_sens.shape)
        # print('the resist fit latent time shape is:\n', fit_t_resist.shape)
        #print('the fit latent time is:\n', fit_t)
        # print('the fit shape is:\n', fit_t.shape)

        return pos, fit_t
    
    def assign_latenttime_v4(self):
        tpoints_sens = self.t1  # t1（0.75）applied to "sensitive" cell
        tpoints_resist = self.t2  # t2（0.25） applied to "resistant" cell

        z_dnn_sens = self.net_z_v2(tpoints_sens)
        z_dnn_resist = self.net_z_v2(tpoints_resist)
        z_obs = self.TGs_expr

        # obtain group information
        groups = self.adata.obs["group"]

        # assign latent time for sensitive cell
        sensitive_mask = (groups == "sensitive").values
        # print('the sensitive mask is:\n', sensitive_mask)
        if sensitive_mask.sum() > 0:
            z_obs_sens = z_obs[sensitive_mask]
            loss_cell_to_t_sens = torch.sum((z_dnn_sens.unsqueeze(1) - z_obs_sens.unsqueeze(0)) ** 2, dim=2)
            pos_sens = torch.argmin(loss_cell_to_t_sens, dim=0)
            pos_sens = torch.clamp(pos_sens, 0, tpoints_sens.size(0) - 1)  # 确保索引不越界
            fit_t_sens = tpoints_sens[pos_sens]
            fit_t_sens = fit_t_sens.flatten()[:, None].squeeze()
        else:
            fit_t_sens = None
            pos_sens = None

        # assign latent time for resistant cell
        resist_mask = (groups == "resistant").values
        # print('the resist mask is:\n', resist_mask)
        if resist_mask.sum() > 0:
            z_obs_resist = z_obs[resist_mask]
            loss_cell_to_t_resist = torch.sum((z_dnn_resist.unsqueeze(1) - z_obs_resist.unsqueeze(0)) ** 2, dim=2)
            pos_resist = torch.argmin(loss_cell_to_t_resist, dim=0)
            pos_resist = torch.clamp(pos_resist, 0, tpoints_resist.size(0) - 1)  # 确保索引不越界
            fit_t_resist = tpoints_resist[pos_resist]
            fit_t_resist = fit_t_resist.flatten()[:, None].squeeze()
        else:
            fit_t_resist = None
            pos_resist = None

        fit_t = torch.zeros(z_obs.shape[0], dtype=tpoints_sens.dtype)
        pos = torch.zeros(z_obs.shape[0], dtype=torch.long)

        if fit_t_sens is not None:
            fit_t[sensitive_mask] = fit_t_sens
            pos[sensitive_mask] = pos_sens
        if fit_t_resist is not None:
            fit_t[resist_mask] = fit_t_resist
            pos[resist_mask] = pos_resist
            
        # print('the sensitive fit latent time shape is:\n', fit_t_sens.shape)
        # print('the resist fit latent time shape is:\n', fit_t_resist.shape)
        print('the fit latent time is:\n', fit_t)
        # print('the fit shape is:\n', fit_t.shape)

        return pos, fit_t

    def calculate_initial_y0(self):
        # calculate initial y0
        V1 = self.V1
        K1 = self.K1
        iroot = self.iroot
        TFLR_allscore = self.TFLR_allscore
        TFs_expr = self.TFs_expr
        # calculate initial y0
        x0 = TFLR_allscore[iroot,:,:]
        Y0 = TFs_expr[iroot,:]
        zero_y = torch.zeros(self.N_TFs, self.N_LRs).float().to(device)
        V1_ = torch.where(x0 > 0, V1, zero_y)
        K1_ = torch.where(x0 > 0, K1, zero_y)
        y0 = torch.sum((V1_ * x0) / ((K1_ + x0) + (1e-12)),dim=1) * Y0
        return y0

    def hill_fun(self, y0, cell_i, t_i):  # trapezoidal rule approximation
        V1 = self.V1
        K1 = self.K1
        TFLR_allscore = self.TFLR_allscore
        TFs_expr = self.TFs_expr
        x_i = TFLR_allscore[int(cell_i), :, :]
        Y_i = TFs_expr[int(cell_i), :]
        zero_y = torch.zeros(self.N_TFs, self.N_LRs)
        V1_ = torch.where(x_i > 0, V1, zero_y)
        K1_ = torch.where(x_i > 0, K1, zero_y)
        tmp1 = torch.sum((V1_ * x_i) / ((K1_ + x_i) + (1e-12)), dim=1) * Y_i
        tmp2 = tmp1 * torch.exp(t_i)
        y_i = (((y0 + tmp2)*t_i)/2 + y0) * torch.exp(-t_i)
        return y_i

    def solve_ym(self, fit_t):
        y0_ = self.calculate_initial_y0()
        N_cell = self.N_cell
        N_TFs = self.N_TFs
        y_ode = torch.zeros((N_cell,N_TFs)).to(device)
        for i in range(N_cell):
            t_i = fit_t[i]
            if t_i.item() == 0:
                y_ode[i] = y0_
            else:
                y_ode[i] = self.hill_fun(y0_,i,t_i)
        return y_ode

    def net_f2(self, isSen2Res):
        N_cell = self.N_cell
        V2 = self.V2
        K2 = self.K2
        regulate = self.regulate
        N_TGs = self.N_TGs
        N_TFs = self.N_TFs
        z_dnn, dz_dt = self.net_z()
        if isSen2Res == True:
            fit_t_pos, fit_t = self.assign_latenttime_v3()
        else:
            fit_t_pos, fit_t = self.assign_latenttime_v4()

        # calculate ym
        y_ode = self.solve_ym(fit_t)

        zero_z = torch.zeros(N_TGs, N_TFs)
        V2_ = torch.where(regulate == 1, V2, zero_z)
        K2_ = torch.where(regulate == 1, K2, zero_z)
        tmp1 = V2_.unsqueeze(0) * y_ode.unsqueeze(1)
        tmp2 = (K2_.unsqueeze(0) + y_ode.unsqueeze(1)) + (1e-12)
        tmp3 = torch.sum(tmp1 / tmp2, dim=2)

        z_pred_exp = torch.zeros((N_cell, N_TGs)).to(device)
        dz_dt_pred = torch.zeros((N_cell, N_TGs)).to(device)
        for i in range(N_cell):
            z_pred_exp[i, :] = z_dnn[fit_t_pos[i]]
            dz_dt_pred[i, :] = dz_dt[fit_t_pos[i]]

        dz_dt_ode = tmp3 - z_pred_exp
        f = dz_dt_pred - dz_dt_ode

        return z_pred_exp, f

    def pre_velo(self,y_ode):
        N_cell = self.N_cell
        N_TGs = self.N_TGs
        # y_ode = self.solve_ym(t)
        pre_velo = torch.zeros((N_cell, N_TGs)).to(device)
        for i in range(N_cell):
            y_i = y_ode[i, :]
            ym_ = self.regulate * y_i
            tmp1 = self.V2 * ym_
            tmp2 = (self.K2 + ym_) + (1e-12)
            tmp3 = torch.sum(tmp1 / tmp2, dim=1)
            dz_dt = tmp3 - self.TGs_expr[i, :]
            pre_velo[i, :] = dz_dt
        return pre_velo

    def train(self, nIter, isSen2Res):
        print('Training SingleVelocity model...')
        self.dnn.train()
        loss_adam = []
        iteration_adam = []
        a = 0
        for epoch in range(nIter):
            z_pred, f_pred = self.net_f2(isSen2Res)
            loss1 = torch.mean((self.TGs_expr - z_pred) ** 2)
            loss2 = torch.mean(f_pred ** 2)  
            
            theta = torch.cat((self.V1.flatten(), self.K1.flatten(), self.V2.flatten(), self.K2.flatten()))
            loss3 = torch.norm(theta, p=1)
            loss4 = torch.norm(theta, p=2)
            # loss3 = torch.norm(self.V1, p=1) + torch.norm(self.K1, p=1) + torch.norm(self.V2, p=1) + torch.norm(self.K2, p=1)
            # loss4 = torch.norm(self.V1, p=2) + torch.norm(self.K1, p=2) + torch.norm(self.V2, p=2) + torch.norm(self.K2, p=2)

            #loss = 0.1 * torch.mean((self.TGs_expr - z_pred) ** 2) + self.Lambda * torch.mean(f_pred ** 2) 
            loss = self.Weight[0] * loss1 + self.Weight[1] * loss2 + self.Weight[2] * loss3 + self.Weight[3] * loss4

            # Backward and optimize
            self.optimizer_Adam.zero_grad()
            loss.backward()
            self.optimizer_Adam.step()
            iteration_adam.append(a)
            a += 1
            loss_adam.append(loss.item())
            # print('It: %d, Loss: %.3e' % (epoch, loss.item()))

            print('loss1: %.3e, loss2: %.3e, loss3: %.3e, loss4: %.3e'%(self.Weight[0] * loss1.item(), 
                                                                        self.Weight[1] * loss2.item(),
                                                                        self.Weight[2] * loss3.item(), 
                                                                        self.Weight[3] * loss4.item()))
            print('It: %d, Loss: %.3e' %(epoch, loss.item()))

        return iteration_adam, loss_adam


def Calculate_Velocity_with_GemOut(adata_copy, model):
    N_TGs = model.N_TGs
    N_TFs = model.N_TFs
    N_cell = model.N_cell
    regulate = model.regulate
    TFs_expr = model.TFs_expr
    TGs_expr = model.TGs_expr
    TGs_pred = model.net_f2(isSen2Res=True)[0]
    V1 = model.V1.detach()
    K1 = model.K1.detach()
    V2 = model.V2.detach()
    K2 = model.K2.detach()

    fit_t = model.assign_latenttime_v3()[1]
    y_ode = model.solve_ym(fit_t)

    adata_copy.obs['latent_time'] = fit_t.detach().numpy()
    print('the latent time is:', fit_t)

    velo_raw = torch.zeros((N_cell, N_TGs)).to(device)
    print('the velo raw shape is:', velo_raw.shape)

    velo_raw_TG = torch.zeros((N_cell, N_TGs)).to(device)
    velo_raw_TF = torch.zeros((N_cell, N_TFs)).to(device)
    print(velo_raw_TG.shape)
    print(velo_raw_TF.shape)

    # calculate Outflow velocity
    for i in range(N_cell):
        y_i = y_ode[i, :]
        ym_ = regulate * y_i
        tmp1 = V2 * ym_
        tmp2 = (K2 + ym_) + (1e-6)
        tmp3 = torch.sum(tmp1 / tmp2, dim=1)
        dz_dt = tmp3 - TGs_expr[i, :]
        velo_raw_TG[i, :] = dz_dt

        # calculate GEM velocity
    for i in range(N_cell):
        y_i = y_ode[i, :]
        ym_ = regulate * y_i
        tmp1 = V2 * ym_
        tmp2 = (K2 + ym_) + (1e-6)
        tmp3 = torch.sum(tmp1 / tmp2, dim=0)
        dz_dt = tmp3 - TFs_expr[i, :]
        velo_raw_TF[i, :] = dz_dt

    velo_raw_Gem = torch.zeros((adata_copy.shape))
    velo_raw_Outflow = torch.zeros((adata_copy.shape))

    GEMs_mask = adata_copy.var['GEMs'].astype(bool)
    GEMs_index = GEMs_mask[GEMs_mask].index
    GEMs_index = [adata_copy.var_names.get_loc(ind) for ind in GEMs_index]  # Convert to integer indices

    Outflows_mask = adata_copy.var['Outflows'].astype(bool)
    Outflows_index = Outflows_mask[Outflows_mask].index
    Outflows_index = [adata_copy.var_names.get_loc(ind) for ind in Outflows_index]  # Convert to integer indices

    print(f"Outflows_index: {Outflows_index}")
    print(f"GEMs_index: {GEMs_index}")

    for i, ind in enumerate(GEMs_index):
        velo_raw_Gem[:, ind] = velo_raw_TF[:, i]

    for i, ind in enumerate(Outflows_index):
        velo_raw_Outflow[:, ind] = velo_raw_TG[:, i]

    print(f"velo_raw_Gem: {velo_raw_Gem}")
    print(f"velo_raw_Gem.shape: {velo_raw_Gem.shape}")

    velo_raw_GemOut = velo_raw_Gem + velo_raw_Outflow

    print(f"velo_raw_GemOut: {velo_raw_GemOut}")
    print(f"velo_raw_GemOut.shape: {velo_raw_GemOut.shape}")

    adata_copy.layers['velo_Gem'] = velo_raw_Gem.detach().numpy()
    adata_copy.layers['velo_GemOut'] = velo_raw_GemOut.detach().numpy()

    return adata_copy


def plotLossAdam(loss_adam, save_path, fig_name):
    torch.save(loss_adam, save_path + "Loss_adam.pt")
    plt.figure(figsize=(10, 4))
    plt.plot(loss_adam)
    # plt.title('Loss of Adam at time %s'%(timepoint+1))
    plt.xlabel('iteration')
    plt.ylabel('loss of Adam')
    plt.savefig(save_path + fig_name, bbox_inches='tight', dpi=600)
    plt.close()



if __name__ == "__main__":
    #Set global random seed
    SEED = 1167
    set_seed(SEED)
    # Load data file
    adata_raw = sc.read('../../data/scPDAC_Malignant_cells_final.h5ad')
    adata_raw.obs['group'] = adata_raw.obs['condition']

    gene_expression_matrix = adata_raw.X.toarray() if isinstance(adata_raw.X, (scipy.sparse.csr_matrix, scipy.sparse.csc_matrix)) else adata_raw.X
    gene_names = adata_raw.var.index  # gene names
    cell_names = adata_raw.obs.index  # cell names
    df_gene_expression = pd.DataFrame(gene_expression_matrix, index=cell_names, columns=gene_names)
    # df_gene_expression.to_csv("../CCCInputData/Dutcal/count_raw_file.csv")

    # data preparation
    df_count = pd.read_csv('../../resluts/SigflowVelo/adata/intermediate_results/count_file.csv',index_col=0)
    df_imput = pd.read_csv('../../resluts/SigflowVelo/adata/intermediate_results/scImputeOutput/scimpute_count.csv',index_col=0)
    InGem_link = pd.read_csv('../../resluts/SigflowVelo/adata/intermediate_results/Inflow-Gem_link_flie.csv', index_col=0)
    GemOut_link = pd.read_csv('../../resluts/SigflowVelo/adata/intermediate_results/Gem-Outflow_link_flie.csv', index_col=0)

    print(f"df_count: {df_count}")
    print(f"GemOut_link: {GemOut_link}")

    # creat AnnData object
    adata = sc.AnnData(X=df_count.values.astype(np.float64))
    adata.obs_names = df_count.index
    adata.var_names = df_count.columns
    adata.obs['celltype'] = adata_raw.obs["celltype"].copy()
    adata.obs['group'] = adata_raw.obs["group"].copy()
    adata.obsm['X_umap'] = adata_raw.obsm['X_umap'].copy()
    adata.layers['Imputate'] = df_imput.values
    print(f"adata: {adata}")

    Inflows = list(np.unique(InGem_link['Inflow'].values))
    GEMs = list(np.unique(GemOut_link['GEM'].values))
    Outflows = list(np.unique(GemOut_link['Outflow'].values))

    ccc_factors = np.unique(np.hstack((Inflows, GEMs, Outflows)))
    print(f"ccc_factors: {ccc_factors}")

    # Annotate corresponding hierarchy in var with 1
    n_gene = adata.shape[1]
    adata.var['Inflows'] = np.full(n_gene, False, dtype=bool).astype(int)
    adata.var['GEMs'] = np.full(n_gene, False, dtype=bool).astype(int)
    adata.var['Outflows'] = np.full(n_gene, False, dtype=bool).astype(int)

    for gene in list(adata.var_names):
        if gene in Inflows:
            adata.var['Inflows'][gene] = 1
        if gene in GEMs:
            adata.var['GEMs'][gene] = 1
        if gene in Outflows:
            adata.var['Outflows'][gene] = 1

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

    # add TFLR score
    InflowGem_score_flie = "../../resluts/SigflowVelo/adata/intermediate_results/Inflow-Gem_score_flie.csv"
    InflowGem_score_matrix = pd.read_csv(InflowGem_score_flie, header=None, skiprows=1, usecols=range(1, 24)).values

    n_dim2, n_dim3 = InflowGem_score_matrix.shape
    InflowGem_allscore = np.tile(InflowGem_score_matrix, (len(adata.obs_names), 1, 1))

    adata.obsm['InflowGem_signaling_score'] = InflowGem_allscore

    # Normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print(f"the structure of adata: {adata}")
    adata.write_h5ad("../../resluts/SigflowVelo/adata/intermediate_results/flowsig_Malignant_cells_CCCInput.h5ad")

    #model train
    adata = sc.read('../../resluts/SigflowVelo/adata/intermediate_results/flowsig_Malignant_cells_CCCInput.h5ad')
    # Prepare Data
    GEMs_expr = torch.tensor(adata.layers['Imputate'][:, adata.var['GEMs'].astype(bool)])
    Outflows_expr = torch.tensor(adata.layers['Imputate'][:, adata.var['Outflows'].astype(bool)])

    GemOut_regulate = torch.tensor(adata.varm['GemOut_regulate'])
    nonzero_idx = torch.nonzero(torch.sum(GemOut_regulate, dim=1)).squeeze()
    GemOut_regulate = GemOut_regulate[nonzero_idx].float()  # torch.Size([539, 114])

    InGem_allscore = torch.tensor(adata.obsm['InflowGem_signaling_score'])

    # Select a cell as root cell
    adata = Function.root_cell(adata, 'UMAP')
    iroot = torch.tensor(adata.uns['iroot'])
    print('the root cell:', adata.uns['iroot'])

    N_Outflows = Outflows_expr.shape[1]
    print(f"the N_Outflows: {N_Outflows}")

    # layers represents the configuration of a neural network, defining the number of neurons in each layer.
    # This is a list indicating the number of neurons in each hidden layer of the neural network.
    # For example, hidden_dims = [64, 32] means the neural network has two hidden layers,
    # with 64 neurons in the first hidden layer and 32 neurons in the second hidden layer.
    hidden_dims = [32, 32]
    layers = hidden_dims
    layers.insert(0, N_Outflows + 1)
    layers.append(N_Outflows)
    # layers = [26, 64, 32, 25] The input layer has 26 nodes, the hidden layer has hidden_dims nodes, and the output layer has 25 nodes.
    print(f"the layers: {layers}")

    data = [Outflows_expr, GEMs_expr, InGem_allscore, GemOut_regulate, iroot, layers]

    TGs_expr, TFs_expr, TFLR_allscore, TGTF_regulate, iroot, layers = data

    # Set hyperparameters
    lr = 0.001
    Weight = [0.1, 0.1, 0.0001, 0.001]
    Lambda = 0.1

    # Initialize model
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model_Sen2Res = SingleVelocity_v2(TGs_expr, TFs_expr, TFLR_allscore, TGTF_regulate, iroot, layers, lr, Weight,
                                      device, adata)

    # Save model and data
    torch.save(model_Sen2Res,"../../resluts/SigflowVelo/adata/intermediate_results/model_Malignant_cells_Sen2Res_v2.pth")
    torch.save(data,"../../resluts/SigflowVelo/adata/intermediate_results/cccvelo_Malignant_cells_v2.pt")

    # Train model
    model_Sen2Res = torch.load("../../resluts/SigflowVelo/adata/intermediate_results/model_Malignant_cells_Sen2Res_v2.pth",weights_only=False)
    data = torch.load("../../resluts/SigflowVelo/adata/intermediate_results/cccvelo_Malignant_cells_v2.pt",weights_only=False)
    adata = sc.read_h5ad("../../resluts/SigflowVelo/adata/intermediate_results/flowsig_Malignant_cells_CCCInput.h5ad")

    # Unpack data
    TGs_expr, TFs_expr, TFLR_allscore, TGTF_regulate, iroot, layers = data

    # lr = 0.001
    # Weight = [0.1, 0.1, 0.0001, 0.001]
    # Lambda = 0.1

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Initialize model
    model_Sen2Res = SingleVelocity_v2(
        TGs_expr=TGs_expr,
        TFs_expr=TFs_expr,
        TFLR_allscore=TFLR_allscore,
        TGTF_regulate=TGTF_regulate,
        iroot=iroot,
        layers=layers,
        lr=lr,
        Weight=Weight,
        device=device,
        adata=adata
    )

    # Train model
    iteration_adam_Sen2Res, loss_adam_Sen2Res = model_Sen2Res.train(nIter=300, isSen2Res=True)
    # Save model after training
    torch.save(model_Sen2Res,"../../resluts/SigflowVelo/adata/intermediate_results/model_Malignant_cells_Sen2Res_trained_v3.pth")
    plotLossAdam(loss_adam_Sen2Res, '../../resluts/SigflowVelo/adata/intermediate_results/','LossAdam_Sen2Res.png')

    # calculate velocity
    model = torch.load("../../resluts/SigflowVelo/adata/intermediate_results/model_Malignant_cells_Sen2Res_trained_v3.pth",weights_only=False)
    adata_copy = adata.copy()
    adata_Sen2Res = Calculate_Velocity_with_GemOut(adata_copy, model)
    adata_Sen2Res.write_h5ad("../../resluts/SigflowVelo/adata/adata_with_calculated_velocity.h5ad")

    ###model_Malignant_cells_Sen2Res_v2.pth -- Initialized model
    ###cccvelo_Malignant_cells_v2.pt  --  Input data: data=TGs_expr, TFs_expr, TFLR_allscore, TGTF_regulate, iroot, layers
    ###model_Malignant_cells_Sen2Res_trained_v3.pth  --  Trained model
    ###adata_with_calculated_velocity.h5ad  --  AnnData with calculated velocity, velocity stored in adata_copy.layers['velo_GemOut']
