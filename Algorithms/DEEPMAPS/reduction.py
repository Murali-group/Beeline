import torch
import torch.utils.data as data
from torch import nn, optim
from torch.nn import functional as F
import numpy as np


class AE(nn.Module):
    def __init__(self, dim):
        super(AE, self).__init__()
        self.dim = dim
        self.fc1 = nn.Linear(dim, 512)
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, 512)
        self.fc4 = nn.Linear(512, dim)

    def encode(self, x):
        h1 = F.relu(self.fc1(x))
        return F.relu(self.fc2(h1))

    def decode(self, z):
        h3 = F.relu(self.fc3(z))
        return torch.relu(self.fc4(h3))

    def forward(self, x):
        z = self.encode(x.view(-1, self.dim))
        return self.decode(z), z


def reduction(method,gene_cell,device):
    if (method == 'AE'):
        encoded,encoded2=reduction_AE(gene_cell,device)
    elif (method == 'raw'):
        encoded = torch.tensor(gene_cell, dtype=torch.float32).to(device)
        encoded2 = torch.tensor(np.transpose(gene_cell),
                                dtype=torch.float32).to(device)
    return encoded,encoded2

def reduction_AE(gene_cell,device):
    gene = torch.tensor(gene_cell, dtype=torch.float32).to(device)
    if gene_cell.shape[0] < 5000:
        ba = gene_cell.shape[0]
    else:
        ba = 5000
    encoded=train_AE(gene,ba,device)

    if gene_cell.shape[1] < 5000:
        ba = gene_cell.shape[1]
    else:
        ba = 5000
    cell = torch.tensor(np.transpose(gene_cell),
                        dtype=torch.float32).to(device)
    encoded2=train_AE(cell,ba,device)
    return encoded,encoded2


def train_AE(feature,ba,device):
    loader = data.DataLoader(feature, ba)
    model = AE(dim=feature.shape[1]).to(device)
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    loss_func = nn.MSELoss()
    EPOCH_AE = 2000
    for epoch in range(EPOCH_AE):
        embedding1 = []
        for _, batch_x in enumerate(loader)	:
            decoded, encoded = model(batch_x)
            loss = loss_func(batch_x, decoded)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            embedding1.append(encoded)
    print('Epoch :', epoch, '|', 'train_loss:%.12f' % loss.data)
    if feature.shape[0] % ba != 0:
        torch.stack(embedding1[0:int(feature.shape[0]/ba)])
        a = torch.stack(embedding1[0:int(feature.shape[0]/ba)])
        a = a.view(ba*int(feature.shape[0]/ba), 256)
        encoded = torch.cat((a, encoded), 0)
    else:
        encode = torch.stack(embedding1)
        encoded = encode.view(feature.shape[0], 256)
    return encoded