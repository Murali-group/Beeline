import os

import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.optim as optim
from torch.autograd import Variable
from torch.utils.data import DataLoader
from torch.utils.data.dataset import TensorDataset

from src.Model import VAE_EAD

Tensor = torch.cuda.FloatTensor


class deepsem_embed:
    def __init__(self, opt):
        self.opt = opt
        try:
            os.mkdir(opt.save_name)
        except:
            print('dir exist')

    def initalize_A(self, data):
        num_genes = data.shape[1]
        A = np.ones([num_genes, num_genes]) / (num_genes - 1) + (np.random.rand(num_genes * num_genes) * 0.0002).reshape(
            [num_genes, num_genes])
        for i in range(len(A)):
            A[i, i] = 0
        return A

    def init_data(self):
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
        data_values = np.maximum(data_values, -15)
        data_values = np.minimum(data_values, 15)
        data = pd.DataFrame(data_values, index=list(data.obs_names), columns=gene_name)
        num_genes, num_nodes = data.shape[1], data.shape[0]
        TF_mask = np.zeros([num_genes, num_genes])
        for i, item in enumerate(data.columns):
            for j, item2 in enumerate(data.columns):
                if i == j:
                    continue
                TF_mask[i, j] = 1
        feat_train = torch.FloatTensor(data.values)
        train_data = TensorDataset(feat_train, torch.LongTensor(list(range(len(feat_train)))),
                                   torch.FloatTensor(Dropout_Mask))
        dataloader = DataLoader(train_data, batch_size=self.opt.batch_size, shuffle=True, num_workers=1)
        return dataloader, num_nodes, num_genes, data, TF_mask, gene_name, Dropout_Mask

    def train_model(self):
        opt = self.opt
        dataloader, num_nodes, num_genes, data, TFmask2, gene_name, Dropout_Mask = self.init_data()
        adj_A_init = self.initalize_A(data)
        vae = VAE_EAD(adj_A_init, 1, opt.n_hidden, opt.K).float().cuda()
        optimizer = optim.RMSprop(vae.parameters(), lr=opt.lr)
        optimizer2 = optim.RMSprop([vae.adj_A], lr=opt.lr * 0.2)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=opt.lr_step_size, gamma=opt.gamma)
        vae.train()
        for epoch in range(opt.n_epochs):
            loss_all, mse_rec, loss_kl, data_ids, loss_tfs, loss_sparse = [], [], [], [], [], []
            if epoch % (opt.K1 + opt.K2) < opt.K1:
                vae.adj_A.requires_grad = False
            else:
                vae.adj_A.requires_grad = True
            for i, data_batch in enumerate(dataloader, 0):
                optimizer.zero_grad()
                inputs, data_id, dropout_mask = data_batch
                inputs = Variable(inputs.type(Tensor))
                data_ids.append(data_id.cpu().detach().numpy())
                temperature = max(0.95 ** epoch, 0.5)
                loss, loss_rec, loss_gauss, loss_cat, dec, y, hidden = vae(inputs, dropout_mask=dropout_mask.cuda(),
                                                                           temperature=temperature, opt=opt)
                sparse_loss = opt.alpha * torch.mean(torch.abs(vae.adj_A))
                loss = loss + sparse_loss
                loss.backward()
                mse_rec.append(loss_rec.item())
                loss_all.append(loss.item())
                loss_kl.append(loss_gauss.item() + loss_cat.item())
                loss_sparse.append(sparse_loss.item())
                if epoch % (opt.K1 + opt.K2) < opt.K1:
                    optimizer.step()
                else:
                    optimizer2.step()
            print('epoch:', epoch, 'loss',
                  np.mean(loss_all), 'mse_loss:', np.mean(mse_rec), 'kl_loss:', np.mean(loss_kl), 'sparse_loss:',
                  np.mean(loss_sparse))
            scheduler.step()

        data_ids = []
        embeds = []
        for i, data_batch in enumerate(dataloader, 0):
            optimizer.zero_grad()
            inputs, data_id, dropout_mask = data_batch
            inputs = Variable(inputs.type(Tensor))
            loss, loss_rec, loss_gauss, loss_cat, dec, y, hidden2 = vae(inputs,
                        dropout_mask=dropout_mask.cuda(), temperature=temperature,opt=opt)
            data_ids.append(data_id.detach().numpy())
            embeds.append(hidden2.cpu().detach().numpy())
        data_ids = np.hstack(data_ids)
        embeds = np.vstack(embeds)
        data_id_map = np.zeros(len(data_ids))
        for i, item in enumerate(data_ids):
            data_id_map[item] = i
        embeds = embeds[data_id_map.astype(int)]
        adata = sc.AnnData(embeds)
        adata.write_h5ad(opt.save_name + '/embedding.h5ad')
