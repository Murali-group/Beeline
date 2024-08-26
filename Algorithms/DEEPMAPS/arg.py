import subprocess
import argparse
import torch.utils.data as data
from kneed import KneeLocator
from sklearn.decomposition import IncrementalPCA
from torch_geometric.data import Data
from pyHGT.model import *
import time
import resource
import datetime
import numpy as np
import random
import torch
from torch import nn, optim
from torch.nn import functional as F
import os
import sys
from pyHGT.data import *
from warnings import filterwarnings
filterwarnings("ignore")
seed = 0

random.seed(seed)
# np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
np.random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)

# torch.cuda.manual_seed(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False


list_in = r.list_in
epoch = list_in['epoch']
lr = list_in['lr']
result_dir = list_in['result_dir']
cuda = list_in['cuda']
gene_cell = list_in['cell_gene']
n_hid = int(list_in['n_hid'])
n_heads = int(list_in['n_heads'])
gene_name = list_in['gene_name']
cell_name = list_in['cell_name']
data_type = list_in['data_type']
if (data_type == 'CITE'):
    n_layers = 4 # type=int, default=4, help='Number of GNN layers'
    sample_depth = 4 # type=int, default=4, help='How many numbers to sample the graph'
    sample_width = 8 # type=int, default=8, help='How many nodes to be sampled per layer per type'
    n_batch = 64 # type=int, default=64, help='Number of batch (sampled graphs) for each epoch'
    batch_size = 64 # type=int, default=64, help='Number of output nodes for training'
    #n_hid = 104
    #n_heads = 13
    #lr = 0.1
    #epoch = 50
if (data_type == 'scRNA_scATAC'):
    n_layers = 2 # type=int, default=4, help='Number of GNN layers'
    sample_depth = 4 # type=int, default=4, help='How many numbers to sample the graph'
    sample_width = 8 # type=int, default=8, help='How many nodes to be sampled per layer per type'
    n_batch = 50 # type=int, default=64, help='Number of batch (sampled graphs) for each epoch'
    batch_size = 110 # type=int, default=64, help='Number of output nodes for training'
    #n_hid = 128
    #n_heads = 16
    #lr = 0.1
    #epoch =100
if (data_type == 'multipleRNA'):
    n_layers = 2 # type=int, default=4, help='Number of GNN layers'
    dropout = 0.0# type=float, default=0, help='Dropout ratio'
    sample_depth = 5 # type=int, default=4, help='How many numbers to sample the graph'
    sample_width = 11 # type=int, default=8, help='How many nodes to be sampled per layer per type'
    n_batch = 80 # type=int, default=64, help='Number of batch (sampled graphs) for each epoch'
    batch_size = 60 # type=int, default=64, help='Number of output nodes for training'
    #n_hid = 104
    #n_heads = 13
    #lr = 0.1
    #epoch = 100

parser = argparse.ArgumentParser(description='Training GNN on gene cell graph')
parser.add_argument('--epoch', type=int, default=epoch)
# Result
parser.add_argument('--result_dir', type=str, default = result_dir,
                    help='The address for storing the models and optimization results.')

#parser.add_argument('--input_dir', type=str, default='default.txt',
#                    help='The address for storing the models and optimization results.')
#parser.add_argument('--label_dir', type=str, default='default.txt',
#                    help='The address for storing the models and optimization results.')
# Feature extration
parser.add_argument('--reduction', type=str, default='AE',
                    help='the method for feature extraction, pca, raw')

parser.add_argument('--in_dim', type=int, default=256,
                    help='Number of hidden dimension (AE)')

# GAE
parser.add_argument('--n_hid', type=int, default=n_hid,
                    help='Number of hidden dimension')
parser.add_argument('--n_heads', type=int, default=n_heads,
                    help='Number of attention head')
parser.add_argument('--n_layers', type=int, default=n_layers,
                    help='Number of GNN layers')
parser.add_argument('--dropout', type=float, default=0,
                    help='Dropout ratio')
parser.add_argument('--sample_depth', type=int, default=sample_depth,
                    help='How many numbers to sample the graph')
parser.add_argument('--sample_width', type=int, default=sample_width,
                    help='How many nodes to be sampled per layer per type')
parser.add_argument('--lr', type=float, default=lr,
                    help='learning rate')
parser.add_argument('--n_batch', type=int, default=n_batch,
                    help='Number of batch (sampled graphs) for each epoch')
parser.add_argument('--batch_size', type=int, default=batch_size,
                    help='Number of output nodes for training')
parser.add_argument('--layer_type', type=str, default='hgt',
                    help='the layer type for GAE')
parser.add_argument('--loss', type=str, default='kl',
                    help='the loss for GAE')

parser.add_argument('--factor', type=float, default='0.5',
                    help='the attenuation factor')
parser.add_argument('--patience', type=int, default=5,
                    help='patience')
parser.add_argument('--rf', type=float, default='0.0',
                    help='the weights of regularization')
parser.add_argument('--cuda', type=int, default=cuda,
                    help='cuda 0 use GPU0 else cpu ')

parser.add_argument('--rep', type=str, default='T',
                    help='precision truncation')

parser.add_argument('--AEtype', type=int, default=1,
                    help='AEtype:1 embedding node autoencoder 2:HGT node autoencode')

parser.add_argument('--optimizer', type=str, default='adamw',
                    help='optimizer')
args = parser.parse_args()


#file0 = 'n_batch'+str(args.n_batch)+'_batch_size_'+str(args.batch_size)+'sample_depth_' +\
#    str(args.sample_depth)+'_nheads_'+str(args.n_heads)+'_nlayers_'+str(args.n_layers) +\
#    '_sample_width_'+str(args.sample_width)+'_lr_'+str(args.lr)+'_n_hid_'+str(args.n_hid) +\
#    '_redution_'+str(args.reduction)+'_rf_'+str(args.rf)+'_factor_'+str(args.factor)+'_pacience_'+str(args.patience)\
#    + '_layertype_'+str(args.layer_type)+'_loss_'+str(args.loss) + \
#    '_optimizer_'+str(args.optimizer)+'_dropout_'+str(args.dropout)
file0='102_n_hid_'+str(args.n_hid)+'_nheads_'+str(args.n_heads)+'_nlayers_'+str(args.n_layers)+'_lr_'+str(args.lr)+'n_batch'+str(args.n_batch)+'batch_size'+str(args.batch_size)


print(file0)
gene_dir = args.result_dir+'/gene/'
cell_dir = args.result_dir+'/cell/'
model_dir = args.result_dir+'/model/'
att_dir = args.result_dir+'/att/'


def load_data(path, sep, col_name, row_name):
    f = open(path, 'r')
    sourceInLine = f.readlines()
    f.close()
    gene_cell = []
    for line in sourceInLine:
        temp1 = line.strip('\n')
        temp2 = temp1.split(sep)
        gene_cell.append(temp2)
    if col_name == True:
        cell_name = gene_cell[0]
        del gene_cell[0]
    else:
        cell_name = ''

    if row_name == True:
        gene_cell = np.array(gene_cell)
        gene_name = gene_cell[:, 0]
        gene_cell = gene_cell[:, 1:gene_cell.shape[1]+1]
    else:
        gene_name = ''
    gene_cell = np.array(gene_cell)
    print("The number of gene is {}, The number of cell is {}".format(
        gene_cell.shape[0], gene_cell.shape[1]))
    return(gene_cell, gene_name, cell_name)


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
        return h1

    def decode(self, z):
        h3 = F.relu(self.fc3(z))
        return torch.relu(self.fc4(h3))
        # return torch.relu(self.fc4(z))

    def forward(self, x):
        z = self.encode(x.view(-1, self.dim))
        return self.decode(z), z


def zeroMean(dataMat):
    meanVal = np.mean(dataMat, axis=0)
    newData = dataMat-meanVal
    return newData, meanVal


def percentage2n(eigVals, percentage):
    sortArray = np.sort(eigVals)
    sortArray = sortArray[-1::-1]
    arraySum = sum(sortArray)
    tmpSum = 0
    num = 0
    for i in sortArray:
        tmpSum += i
        num += 1
        if tmpSum >= arraySum*percentage:
            return num


def pca_1(dataMat, percentage=0.99):
    newData, meanVal = zeroMean(dataMat)
    covMat = np.cov(newData, rowvar=0)
    eigVals, eigVects = np.linalg.eig(np.mat(covMat))
    n = percentage2n(eigVals, percentage)
    eigValIndice = np.argsort(eigVals)
    n_eigValIndice = eigValIndice[-1:-(n+1):-1]
    n_eigVect = eigVects[:, n_eigValIndice]
    lowDDataMat = newData*n_eigVect
    reconMat = (lowDDataMat*n_eigVect.T)+meanVal
    return lowDDataMat, reconMat


def pca(dataMat, n):
    newData, meanVal = zeroMean(dataMat)
    covMat = np.cov(newData, rowvar=0)

    eigVals, eigVects = np.linalg.eig(np.mat(covMat))
    eigValIndice = np.argsort(eigVals)
    n_eigValIndice = eigValIndice[-1:-(n+1):-1]
    n_eigVect = eigVects[:, n_eigValIndice]
    lowDDataMat = newData*n_eigVect
    reconMat = (lowDDataMat*n_eigVect.T)+meanVal
    return lowDDataMat, reconMat


def debuginfoStr(info):

    # print('---'+str(datetime.timedelta(seconds=int(time.time()-start_time)))+'---'+info)
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
    print('Mem consumption: '+str(mem))


def sub_sample(samp_nodes, sampled_depth=args.sample_depth, sampled_number=args.sample_width):
    inp = {'gene': np.concatenate(
        [samp_nodes, graph.years[samp_nodes]]).reshape(2, -1).transpose()}
    '''
        Sample Sub-Graph based on the connection of other nodes with currently sampled nodes
        We maintain budgets for each node type, indexed by <node_id, time>.
        Currently sampled nodes are stored in layer_data.
        After nodes are sampled, we construct the sampled adjacancy matrix.
    '''
    layer_data = defaultdict(  # target_type
        lambda: {}  # {target_id: [ser, time]}
    )

    budget = defaultdict(  # source_type
        lambda: defaultdict(  # source_id
            # [sampled_score, time]
            lambda: [0., 0]
        ))
    new_layer_adj = defaultdict(  # target_type
        lambda: defaultdict(  # source_type
            lambda: defaultdict(  # relation_type
                lambda: []  # [target_id, source_id]
            )))
    '''
        For each node being sampled, we find out all its neighborhood, 
        adding the degree count of these nodes in the budget.
        Note that there exist some nodes that have many neighborhoods
        (such as fields, venues), for those case, we only consider 
    '''
    def add_budget(te, target_id, target_time, layer_data, budget):
        for source_type in te:
            tes = te[source_type]
            for relation_type in tes:
                if relation_type == 'self' or target_id not in tes[relation_type]:
                    continue
                adl = tes[relation_type][target_id]
                if len(adl) < sampled_number:
                    sampled_ids = list(adl.keys())
                else:
                    sampled_ids = np.random.choice(
                        list(adl.keys()), sampled_number, replace=False)
                for source_id in sampled_ids:
                    source_time = adl[source_id]
                    if source_time == None:
                        source_time = target_time
                    if source_id in layer_data[source_type]:
                        continue
                    budget[source_type][source_id][0] += 1. / len(sampled_ids)
                    budget[source_type][source_id][1] = source_time

    '''
        First adding the sampled nodes then updating budget.
    '''
    for _type in inp:
        for _id, _time in inp[_type]:
            layer_data[_type][_id] = [len(layer_data[_type]), _time]
    for _type in inp:
        te = graph.edge_list[_type]
        for _id, _time in inp[_type]:
            add_budget(te, _id, _time, layer_data, budget)
    '''
        We recursively expand the sampled graph by sampled_depth.
        Each time we sample a fixed number of nodes for each budget,
        based on the accumulated degree.
    '''
    for layer in range(sampled_depth):
        sts = list(budget.keys())
        for source_type in sts:
            te = graph.edge_list[source_type]
            keys = np.array(list(budget[source_type].keys()))
            if sampled_number > len(keys):
                '''
                    Directly sample all the nodes
                '''
                sampled_ids = np.arange(len(keys))
            else:
                '''
                    Sample based on accumulated degree
                '''
                score = np.array(list(budget[source_type].values()))[:, 0] ** 2
                score = score / np.sum(score)
                sampled_ids = np.random.choice(
                    len(score), sampled_number, p=score, replace=False)
            sampled_keys = keys[sampled_ids]
            '''
                First adding the sampled nodes then updating budget.
            '''
            for k in sampled_keys:
                layer_data[source_type][k] = [
                    len(layer_data[source_type]), budget[source_type][k][1]]
            for k in sampled_keys:

                add_budget(te, k, budget[source_type]
                           [k][1], layer_data, budget)
                budget[source_type].pop(k)

    '''
        Prepare feature, time and adjacency matrix for the sampled graph
    '''
    feature = {}
    times = {}
    indxs = {}
    texts = []
    for _type in layer_data:
        # print(_type)
        if len(layer_data[_type]) == 0:
            continue
        idxs = np.array(list(layer_data[_type].keys()), dtype=np.int)
        # print(idxs)
        tims = np.array(list(layer_data[_type].values()))[:, 1]
        if _type == 'cell':
            idxs = idxs - gene_cell.shape[0]
        feature[_type] = graph.node_feature[_type][idxs]
        times[_type] = tims
        indxs[_type] = idxs

    edge_list = defaultdict(  # target_type
        lambda: defaultdict(  # source_type
            lambda: defaultdict(  # relation_type
                lambda: []  # [target_id, source_id]
            )))
    for _type in layer_data:
        for _key in layer_data[_type]:
            _ser = layer_data[_type][_key][0]
            edge_list[_type][_type]['self'] += [[_ser, _ser]]
    '''
        Reconstruct sampled adjacancy matrix by checking whether each
        link exist in the original graph
    '''
    for target_type in graph.edge_list:
        te = graph.edge_list[target_type]
        tld = layer_data[target_type]
        for source_type in te:
            tes = te[source_type]
            sld = layer_data[source_type]
            for relation_type in tes:
                tesr = tes[relation_type]
                for target_key in tld:
                    if target_key not in tesr:
                        continue
                    target_ser = tld[target_key][0]
                    for source_key in tesr[target_key]:
                        '''
                            Check whether each link (target_id, source_id) exist in original adjacancy matrix
                        '''
                        if source_key in sld:
                            source_ser = sld[source_key][0]
                            edge_list[target_type][source_type][relation_type] += [
                                [target_ser, source_ser]]
    # print("feature",feature)

    return feature, times, edge_list, indxs, texts


# load data
start_time = time.time()
print('---0:00:00---scRNA starts loading.')
#gene_cell, gene_name, cell_name = load_data(
#    args.input_dir, sep=" ", col_name=True, row_name=True)
gene_cell = gene_cell.astype('float')
# gene_cell=gene_cell[0:1000,:]
# cpu/gpu
cuda = args.cuda  # 'cpu'#-1
if cuda == 0:
    device = torch.device("cuda:" + "0")
    print("cuda>>>")
else:
    device = torch.device("cpu")
print(device)


debuginfoStr('scRNA has been successfully loaded')

#args_reduction ='pca'
# feature extraction

if (args.reduction == 'AE'):
    gene = torch.tensor(gene_cell, dtype=torch.float32).to(device)
    if gene_cell.shape[0] < 5000:
        ba = gene_cell.shape[0]
    else:
        ba = 5000
    loader1 = data.DataLoader(gene, ba)

    EPOCH_AE = 2000
    model = AE(dim=gene.shape[1]).to(device)
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    loss_func = nn.MSELoss()
    for epoch in range(EPOCH_AE):
        embedding1 = []
        for _, batch_x in enumerate(loader1)	:

            decoded, encoded = model(batch_x)
        #encoded1 , decoded1 = Coder2(cell)
            loss = loss_func(batch_x, decoded)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            embedding1.append(encoded)
    print('Epoch :', epoch, '|', 'train_loss:%.12f' % loss.data)
    if gene.shape[0] % ba != 0:
        torch.stack(embedding1[0:int(gene.shape[0]/ba)])
        a = torch.stack(embedding1[0:int(gene.shape[0]/ba)])
        a = a.view(ba*int(gene.shape[0]/ba), 256)
        encoded = torch.cat((a, encoded), 0)

    else:
        encode = torch.stack(embedding1)
        encoded = encode.view(gene.shape[0], 256)

    if gene_cell.shape[1] < 5000:
        ba = gene_cell.shape[1]
    else:
        ba = 5000
    cell = torch.tensor(np.transpose(gene_cell),
                        dtype=torch.float32).to(device)
    loader2 = data.DataLoader(cell, ba)
    model2 = AE(dim=cell.shape[1]).to(device)
    optimizer2 = optim.Adam(model2.parameters(), lr=1e-3)
    for epoch in range(EPOCH_AE):
        embedding1 = []
        for _, batch_x in enumerate(loader2):
            decoded2, encoded2 = model2(batch_x)
            loss = loss_func(batch_x, decoded2)
            optimizer2.zero_grad()
            loss.backward()
            optimizer2.step()
            embedding1.append(encoded2)
    print('Epoch :', epoch, '|', 'train_loss:%.12f' % loss.data)
    if cell.shape[0] % ba != 0:
        torch.stack(embedding1[0:int(cell.shape[0]/ba)])
        a = torch.stack(embedding1[0:int(cell.shape[0]/ba)])
        a = a.view(ba*int(cell.shape[0]/ba), 256)
        encoded2 = torch.cat((a, encoded2), 0)
        # encode.shape
    else:
        encode = torch.stack(embedding1)
        encoded2 = encode.view(cell.shape[0], 256)

if (args.reduction == 'raw'):
    encoded = torch.tensor(gene_cell, dtype=torch.float32).to(device)
    encoded2 = torch.tensor(np.transpose(gene_cell),
                            dtype=torch.float32).to(device)


debuginfoStr('Feature extraction finished')


target_nodes = np.arange(gene_cell.shape[1]+gene_cell.shape[0])
# gene cell
g = np.nonzero(gene_cell)[0]
c = np.nonzero(gene_cell)[1]+gene_cell.shape[0]
edge1 = list(g)
edge2 = list(c)
# print(len(edge1))
# print(len(edge2))
#node_feature = torch.cat( (encoded, encoded2), 0)
#node_feature = (node_feature-torch.mean(node_feature))/torch.std(node_feature)
edge_index = torch.tensor([edge1, edge2], dtype=torch.float)
# x={'gene': torch.tensor(node_feature[0:encoded.shape[0],:], dtype=torch.float),
#    'cell': torch.tensor(node_feature[encoded.shape[0]:(encoded2.shape[0]+encoded.shape[0]),:], dtype=torch.float),
# }
x = {'gene': torch.tensor(encoded, dtype=torch.float),
     'cell': torch.tensor(encoded2, dtype=torch.float),
     }
#x = torch.tensor(np.random.randn(len(target_nodes),64), dtype=torch.float)
#h = Data(edge_index_dict=edge_index, x=x)
edge_index_dict = {('gene', 'g_c', 'cell')
                    : torch.tensor([g, c], dtype=torch.float)}
edge_reltype = {
    ('gene', 'g_c', 'cell'):  torch.tensor([g, c]).shape[1]
}
num_nodes_dict = {
    'gene': gene_cell.shape[0],
    'cell': gene_cell.shape[1]
}
data = Data(edge_index_dict=edge_index_dict,
            edge_reltype=edge_reltype, num_nodes_dict=num_nodes_dict, x=x)
graph = Graph()
edg = graph.edge_list
edge_index_dict = data.edge_index_dict
for key in edge_index_dict:
    # print(key)
    edges = edge_index_dict[key]
    s_type, r_type, t_type = key[0], key[1], key[2]
    elist = edg[t_type][s_type][r_type]
    rlist = edg[s_type][t_type]['rev_' + r_type]
    for s_id, t_id in edges.t().tolist():
        year = 1
        elist[t_id][s_id] = year
        rlist[s_id][t_id] = year
edg = {}
deg = {key: np.zeros(data.num_nodes_dict[key]) for key in data.num_nodes_dict}

for k1 in graph.edge_list:
    if k1 not in edg:
        edg[k1] = {}
    for k2 in graph.edge_list[k1]:
        if k2 not in edg[k1]:
            edg[k1][k2] = {}
        for k3 in graph.edge_list[k1][k2]:
            if k3 not in edg[k1][k2]:
                edg[k1][k2][k3] = {}
            for num1, e1 in enumerate(graph.edge_list[k1][k2][k3]):
                if len(graph.edge_list[k1][k2][k3][e1]) == 0:
                    continue

                edg[k1][k2][k3][num1] = {}
                for num2, e2 in enumerate(graph.edge_list[k1][k2][k3][e1]):
                    edg[k1][k2][k3][num1][num2] = graph.edge_list[k1][k2][k3][e1][e2]
                deg[k1][num1] += len(edg[k1][k2][k3][num1])
            #print(k1, k2, k3, len(edg[k1][k2][k3]))

graph.node_feature['gene'] = data.x['gene']
graph.node_feature['cell'] = data.x['cell']
print(graph.node_feature['gene'].shape)
print(graph.node_feature['cell'].shape)


graph.years = np.ones(len(target_nodes))

np.random.seed(seed)
jobs = []

for batch_id in np.arange(args.n_batch):

    p = sub_sample(np.random.choice(
        np.arange(gene_cell.shape[0]), args.batch_size, replace=False))
    jobs.append(p)
debuginfoStr('Cell Graph constructed and pruned')

if (args.reduction != 'raw'):
    gnn = GNN(conv_name=args.layer_type, in_dim=encoded.shape[1],
              n_hid=args.n_hid, n_heads=args.n_heads, n_layers=args.n_layers, dropout=args.dropout,
              num_types=2, num_relations=2, use_RTE=False).to(device)
else:
    gnn = GNN_from_raw(conv_name=args.layer_type, in_dim=[encoded.shape[1], encoded2.shape[1]],
                       n_hid=args.n_hid, n_heads=args.n_heads, n_layers=args.n_layers, dropout=args.dropout,
                       num_types=2, num_relations=2, use_RTE=False,
                       AEtype=args.AEtype).to(device)


#args.optimizer =  'adamw'
if args.optimizer == 'adamw':
    optimizer = torch.optim.AdamW(gnn.parameters(), lr=args.lr)
elif args.optimizer == 'adam':
    optimizer = torch.optim.Adam(gnn.parameters(), lr=args.lr)
elif args.optimizer == 'sgd':
    optimizer = torch.optim.SGD(gnn.parameters(), lr=args.lr)
elif args.optimizer == 'adagrad':
    optimizer = torch.optim.Adagrad(gnn.parameters(), lr=args.lr)
# gnn.double()

#model, optimizer = amp.initialize(gnn, optimizer, opt_level="O1")
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
    optimizer, 'min', factor=args.factor, patience=args.patience, verbose=True)
#scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.9)

gnn.train()
k = []
for epoch in np.arange(args.epoch):
    L = 0
    for job in jobs:
        # print(job)
        feature = job[0]
        time = job[1]
        edge_list = job[2]
        indxs = job[3]
        node_dict = {}
        node_feature = []
        node_type = []
        node_time = []
        edge_index = []
        edge_type = []
        edge_time = []

        node_num = 0
        types = graph.get_types()
        for t in types:
            #print("t in types "+str(t)+"\n")
            node_dict[t] = [node_num, len(node_dict)]
            node_num += len(feature[t])
            if args.reduction == 'raw':
                node_feature.append([])

        # for t_i in range(len(types)):
        for t in types:
            t_i = node_dict[t][1]
            #print("feature t:\n")
            #print("t_i="+str(t_i)+" t="+str(t)+"\n")
            # print(feature[t].shape)
            if args.reduction != 'raw':
                node_feature += list(feature[t])
            else:
                node_feature[t_i] = torch.tensor(
                    feature[t], dtype=torch.float32).to(device)

            node_time += list(time[t])
            node_type += [node_dict[t][1] for _ in range(len(feature[t]))]

        # print("node_type:\n")
        # print(node_type)
        edge_dict = {e[2]: i for i, e in enumerate(graph.get_meta_graph())}
        edge_dict['self'] = len(edge_dict)

        for target_type in edge_list:
            for source_type in edge_list[target_type]:
                for relation_type in edge_list[target_type][source_type]:
                    for ii, (ti, si) in enumerate(edge_list[target_type][source_type][relation_type]):
                        tid, sid = ti + \
                            node_dict[target_type][0], si + \
                            node_dict[source_type][0]
                        edge_index += [[sid, tid]]
                        edge_type += [edge_dict[relation_type]]

                        # Our time ranges from 1900 - 2020, largest span is 120.

                        edge_time += [node_time[tid] - node_time[sid] + 120]

        if (args.reduction != 'raw'):
            node_feature = torch.stack(node_feature)
            node_feature = torch.tensor(node_feature, dtype=torch.float32)
            node_feature = node_feature.to(device)

        #node_feature = torch.trunc(node_feature*10000)/10000
        node_type = torch.LongTensor(node_type)
        edge_time = torch.LongTensor(edge_time)
        edge_index = torch.LongTensor(edge_index).t()
        edge_type = torch.LongTensor(edge_type)
        if (args.reduction == 'raw'):
            node_rep, node_decoded_embedding = gnn.forward(node_feature, node_type.to(device),
                                                           edge_time.to(
                                                               device),
                                                           edge_index.to(
                                                               device),
                                                           edge_type.to(device))
        else:
            node_rep = gnn.forward(node_feature, node_type.to(device),
                                   edge_time.to(device),
                                   edge_index.to(device),
                                   edge_type.to(device))

        if args.rep == 'T':
            node_rep = torch.trunc(node_rep*10000000000)/10000000000
            if args.reduction == 'raw':
                for t in types:
                    t_i = node_dict[t][1]
                    # print("t_i="+str(t_i))
                    node_decoded_embedding[t_i] = torch.trunc(
                        node_decoded_embedding[t_i]*10000000000)/10000000000

        # print(node_rep)
        # print(abc)
        gene_matrix = node_rep[node_type == 0, ]
        cell_matrix = node_rep[node_type == 1, ]
        # print("gene_matrix shape \n") #(size, 64)
        # print(gene_matrix.shape)
        #print("cell_matrix shape \n") (size, 64)
        # print(cell_matrix.shape)
        # args_loss='kl'
        regularization_loss = 0
        for param in gnn.parameters():
            regularization_loss += torch.sum(torch.pow(param, 2))
        if (args.loss == "kl"):
            decoder = torch.mm(gene_matrix, cell_matrix.t())
            adj = gene_cell[indxs['gene'], ]
            adj = adj[:, indxs['cell']]
            adj = torch.tensor(adj, dtype=torch.float32).to(device)
            if args.reduction == 'raw':
                if epoch % 2 == 0:
                    loss = F.kl_div(decoder.softmax(
                        dim=-1).log(), adj.softmax(dim=-1), reduction='sum')+args.rf*regularization_loss
                else:
                    loss = nn.MSELoss()(
                        node_feature[0], node_decoded_embedding[0])+args.rf*regularization_loss
                    # print("loss=\n")
                    # print(loss)
                    for t_i in range(1, len(types)):
                        loss += nn.MSELoss()(node_feature[t_i],
                                             node_decoded_embedding[t_i])
            else:
                loss = F.kl_div(decoder.softmax(dim=-1).log(),
                                adj.softmax(dim=-1), reduction='sum')

            # print("loss1",loss)
            #loss = torch.trunc(loss*10000)/100
            # loss.requires_grad_(True)
            # print("loss2",loss)
        if (args.loss == "cross"):
            from torch_geometric.utils import remove_self_loops
            from torch_geometric.utils import add_self_loops

            EPS = 1e-15
            value = (node_rep[edge_index[0]] *
                     node_rep[edge_index[1]]).sum(dim=1)
            torch.sigmoid(value)
            pos_loss = -torch.log(
                torch.sigmoid(value) + EPS).mean()
            pos_edge_index = edge_index

            pos_edge_index, _ = remove_self_loops(pos_edge_index)
            pos_edge_index, _ = add_self_loops(pos_edge_index)
            neg_edge_index = negative_sampling(
                pos_edge_index, node_rep.size(0))
            value = (node_rep[neg_edge_index[0]] *
                     node_rep[neg_edge_index[1]]).sum(dim=1)
            neg_loss = -torch.log(1 - torch.sigmoid(value) + EPS).mean()

            loss = neg_loss + pos_loss
            # loss=torch.trunc(loss*1000000000)/1000000000

        L += loss.item()
        # print(L)
        #loss.backward( retain_graph = True)
        # loss.requires_grad()
        optimizer.zero_grad()
        # with amp.scale_loss(loss, optimizer) as scaled_loss:
        #   scaled_loss.backward()
        loss.backward()
        # print(loss.grad)
        # torch.nn.utils.clip_grad_norm_(gnn.parameters(),1)
        optimizer.step()
    # print(node_rep)
    scheduler.step(L/(int(gene_cell.shape[0])))
    # k.append(L/(int(gene_cell.shape[0]/10)))
    # print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
    print('Epoch :', epoch+1, '|', 'train_loss:%.12f' %
          (L/(int(gene_cell.shape[0]))/args.n_batch))


state = {'model': gnn.state_dict(), 'optimizer': scheduler.state_dict(),
         'epoch': epoch}
torch.save(state, model_dir+file0)


debuginfoStr('Graph Autoencoder training finished')

debuginfoStr('load training model')

state = torch.load(model_dir+file0, map_location=lambda storage, loc: storage)
device = torch.device("cpu")

if (args.reduction != 'raw'):
    gnn = GNN(conv_name=args.layer_type, in_dim=encoded.shape[1], n_hid=args.n_hid, n_heads=args.n_heads, n_layers=args.n_layers, dropout=args.dropout,
              num_types=2, num_relations=2, use_RTE=False).to(device)
else:
    gnn = GNN_from_raw(conv_name=args.layer_type, in_dim=[encoded.shape[1], encoded2.shape[1]], n_hid=args.n_hid, n_heads=args.n_heads, n_layers=args.n_layers, dropout=args.dropout,
                       num_types=2, num_relations=2, use_RTE=False,
                       AEtype=args.AEtype).to(device)


# model.eval()
#if (gene_cell.shape[0] > 10000):
#    ba = 5000
#else:
#    ba = gene_cell.shape[0]
if (gene_cell.shape[1]>10000):

    if (gene_cell.shape[0]>10000):
        ba = 500
    else:
        ba = gene_cell.shape[0]
else:
    if (gene_cell.shape[0]>10000):
        ba = 5000
    else:
        ba = gene_cell.shape[0]

#optimizer = torch.optim.Adagrad(model.parameters(),lr=args.lr)
scheduler.load_state_dict(state['optimizer'])
gnn.load_state_dict(state['model'])
g_embedding = []
gene_name = []
cell_name = []
attention = []
with torch.no_grad():
    for i in range(0, gene_cell.shape[0], ba):
        adj = gene_cell[i:(i+ba), :]  # number of gene x all cells
        g = np.nonzero(adj)[0]
        c = np.nonzero(adj)[1]+adj.shape[0]
        edge1 = list(g)
        edge2 = list(c)
        edge_index = torch.tensor([edge1, edge2], dtype=torch.long)
        x = {'gene': torch.tensor(encoded[i:(ba+i), :], dtype=torch.float),
             'cell': torch.tensor(encoded2, dtype=torch.float),
             }  # batch of gene all cells
        edge_index_dict = {('gene', 'g_c', 'cell')
                            : torch.tensor([g, c], dtype=torch.long)}
        edge_reltype = {('gene', 'g_c', 'cell')
                         :  torch.tensor([g, c]).shape[1]}
        num_nodes_dict = {'gene': adj.shape[0], 'cell': gene_cell.shape[1]}
        data = Data(edge_index_dict=edge_index_dict,
                    edge_reltype=edge_reltype, num_nodes_dict=num_nodes_dict, x=x)
        a = np.nonzero(adj)[0]
        b = np.nonzero(adj)[1]
        node_type = list(np.zeros(adj.shape[0]))+list(np.ones(adj.shape[1]))
        node_type = torch.LongTensor(node_type)
        edge_index = data['edge_index_dict'][('gene',  'g_c',  'cell')]
        edge_type = list(np.zeros(len(edge_index[1])))
        edge_time = torch.LongTensor(list(np.zeros(len(edge_index[1]))))
        edge_type = torch.LongTensor(edge_type)
        print(len(x['gene']))  # 5000
        print(len(x['cell']))  # 2713
        print(len(node_type))  # 7713
        if args.reduction != 'raw':
            node_rep = gnn.forward((torch.cat((x['gene'], x['cell']), 0)).to(device), node_type.to(device),
                                   edge_time.to(device), edge_index.to(device), edge_type.to(device))
        else:
            node_rep, _ = gnn.forward([x['gene'].to(device), x['cell'].to(device)], node_type.to(device),
                                      edge_time.to(device), edge_index.to(device), edge_type.to(device))

        gene_name = gene_name + list(np.array(edge_index[0]+i))
        cell_name = cell_name + list(np.array(edge_index[1]-adj.shape[0]))
        attention.append(gnn.att)
        gene_matrix = node_rep[node_type == 0, ]
        cell_matrix = node_rep[node_type == 1, ]
        g_embedding.append(gene_matrix)

if gene_cell.shape[0] % ba == 0:
    gene_matrix = np.vstack(g_embedding[0:int(gene_cell.shape[0]/ba)])
    attention = np.vstack(attention[0:int(gene_cell.shape[0]/ba)])
else:
    final_tensor = np.vstack(g_embedding[0:int(gene_cell.shape[0]/ba)])
    gene_matrix = np.concatenate((final_tensor, gene_matrix), 0)
    final_attention = np.vstack(attention[0:int(gene_cell.shape[0]/ba)])
    attention = np.concatenate((final_attention, gnn.att), 0)
cell_matrix = cell_matrix.detach().numpy()
#np.savetxt(gene_dir+file0, gene_matrix, delimiter=' ')
#np.savetxt(cell_dir+file0, cell_matrix, delimiter=' ')
debuginfoStr(' finished')
g = np.nonzero(gene_cell)[0]
c = np.nonzero(gene_cell)[1]+gene_cell.shape[0]
name1 = pd.DataFrame(
    gene_name[0:torch.tensor([g, c]).shape[1]], columns=['gene'])
name2 = pd.DataFrame(
    cell_name[0:torch.tensor([g, c]).shape[1]], columns=['cell'])
df = pd.DataFrame(attention)
df2 = pd.concat([name1, name2, df], axis=1)
attention = df2
#df2.to_csv(att_dir+file0, sep=",", index=True)
