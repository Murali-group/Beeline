import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import scanpy as sc
import torch.optim as optim
from torch.autograd import Variable
from torch.utils.data import DataLoader
from torch.utils.data.dataset import TensorDataset
from src.utils import evaluate, extractEdgesFromMatrix
from src.Model import VAE_EAD


class celltype_GRN_model:
    def __init__(self,opt):
        self.opt = opt
        try:
            os.mkdir(opt.save_name)
        except:
            print('dir exist')
    def initalize_A(self,data):
        num_genes = data.shape[1]
        A = np.ones([num_genes, num_genes]) / (num_genes - 1) + (np.random.rand(num_genes * num_genes) * 0.0002).reshape(
            [num_genes, num_genes])
        for i in range(len(A)):
            A[i, i] = 0
        return A


    def init_data(self,):
        Ground_Truth = pd.read_csv(self.opt.net_file, header=0)
        data = sc.read(self.opt.data_file)
        gene_name = list(data.var_names)
        data_values = data.X
        Dropout_Mask = (data_values != 0).astype(float)
        means = []
        stds = []
        for i in range(data_values.shape[1]):
            tmp = data_values[:, i]
            means.append(tmp[tmp != 0].mean())
            stds.append(tmp[tmp != 0].std())
        means = np.array(means)
        stds = np.array(stds)
        stds[np.isnan(stds)] = 0
        stds[np.isinf(stds)] = 0
        data_values = (data_values - means) / (stds)
        data_values[np.isnan(data_values)] = 0
        data_values[np.isinf(data_values)] = 0
        data_values = np.maximum(data_values, -10)
        data_values = np.minimum(data_values, 10)
        data = pd.DataFrame(data_values, index=list(data.obs_names), columns=gene_name)
        TF = set(Ground_Truth['Gene1'])
        All_gene = set(Ground_Truth['Gene1']) | set(Ground_Truth['Gene2'])
        num_genes, num_nodes = data.shape[1], data.shape[0]
        Evaluate_Mask = np.zeros([num_genes, num_genes])
        TF_mask = np.zeros([num_genes, num_genes])
        for i, item in enumerate(data.columns):
            for j, item2 in enumerate(data.columns):
                if i == j:
                    continue
                if item2 in TF and item in All_gene:
                    Evaluate_Mask[i, j] = 1
                if item2 in TF:
                    TF_mask[i, j] = 1
        feat_train = torch.FloatTensor(data.values)
        train_data = TensorDataset(feat_train, torch.LongTensor(list(range(len(feat_train)))),
                                   torch.FloatTensor(Dropout_Mask))
        dataloader = DataLoader(train_data, batch_size=self.opt.batch_size, shuffle=True, num_workers=1)
        truth_df = pd.DataFrame(np.zeros([num_genes, num_genes]), index=data.columns, columns=data.columns)
        for i in range(Ground_Truth.shape[0]):
            truth_df.loc[Ground_Truth.iloc[i, 1], Ground_Truth.iloc[i, 0]] = 1
        A_truth = truth_df.values
        idx_rec, idx_send = np.where(A_truth)
        truth_edges = set(zip(idx_send, idx_rec))
        return dataloader, Evaluate_Mask, num_nodes, num_genes, data, truth_edges, TF_mask, gene_name,


    def train_model(self):
        dataloader, Evaluate_Mask, num_nodes, num_genes, data, truth_edges, TFmask2, gene_name  = self.init_data()
        adj_A_init  = self.initalize_A(data)
        vae = VAE_EAD(adj_A_init, 1, self.opt.n_hidden, self.opt.K).float().cuda()
        Tensor = torch.cuda.FloatTensor
        optimizer = optim.RMSprop(vae.parameters(), lr=self.opt.lr)
        optimizer2 = optim.RMSprop([vae.adj_A], lr=self.opt.lr * 0.2)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=self.opt.lr_step_size, gamma=self.opt.gamma)
        best_Epr = 0
        vae.train()
        print(vae)
        for epoch in range(self.opt.n_epochs+1):
            loss_all, mse_rec, loss_kl, data_ids, loss_tfs, loss_sparse = [], [], [], [], [], []
            if epoch % (self.opt.K1+self.opt.K2) < self.opt.K1:
                vae.adj_A.requires_grad = False
            else:
                vae.adj_A.requires_grad = True
            for i, data_batch in enumerate(dataloader, 0):
                optimizer.zero_grad()
                inputs, data_id, dropout_mask = data_batch
                inputs = Variable(inputs.type(Tensor))
                data_ids.append(data_id.cpu().detach().numpy())
                temperature = max(0.95 ** epoch, 0.5)
                loss, loss_rec, loss_gauss, loss_cat, dec, y, hidden = vae(inputs,dropout_mask=dropout_mask.cuda(),temperature=temperature,opt=self.opt)
                sparse_loss = self.opt.alpha * torch.mean(torch.abs(vae.adj_A))
                loss = loss + sparse_loss
                loss.backward()
                mse_rec.append(loss_rec.item())
                loss_all.append(loss.item())
                loss_kl.append(loss_gauss.item() + loss_cat.item())
                loss_sparse.append(sparse_loss.item())
                if epoch % (self.opt.K1+self.opt.K2) < self.opt.K1:
                    optimizer.step()
                else:
                    optimizer2.step()
            scheduler.step()
            if epoch % (self.opt.K1+self.opt.K2) >= self.opt.K1:
                Ep, Epr =  evaluate(vae.adj_A.cpu().detach().numpy(), truth_edges, Evaluate_Mask)
                best_Epr = max(Epr, best_Epr)
                print('epoch:', epoch, 'Ep:', Ep, 'Epr:', Epr, 'loss:',
                      np.mean(loss_all), 'mse_loss:', np.mean(mse_rec), 'kl_loss:', np.mean(loss_kl), 'sparse_loss:',np.mean(loss_sparse))
        extractEdgesFromMatrix(vae.adj_A.cpu().detach().numpy(), gene_name,TFmask2).to_csv(
        self.opt.save_name + '/GRN_inference_result.tsv', sep='\t', index=False)

