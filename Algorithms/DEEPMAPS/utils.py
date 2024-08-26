from pyHGT.data import Graph
import numpy as np
import torch
from collections import defaultdict
import resource
import pandas as pd

def debuginfoStr(info):
    print(info)
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024*1024)
    print('Mem consumption (GB): '+str(mem))


def loadGAS(data_path):
    df=pd.read_csv(data_path, sep=" ")
    return df.to_numpy()



def build_graph(gene_cell,encoded,encoded2):
    g_index,c_index = np.nonzero(gene_cell)
    # 加上偏移量作为cell的节点标号
    c_index += gene_cell.shape[0]
    edges=torch.tensor([g_index, c_index], dtype=torch.float)
    
    # 这里是直接对graph.edge_list进行修改了，不是副本
    graph = Graph()
    s_type, r_type, t_type = ('gene', 'g_c', 'cell')
    elist = graph.edge_list[t_type][s_type][r_type]
    rlist = graph.edge_list[s_type][t_type]['rev_' + r_type]
    year = 1
    for s_id, t_id in edges.t().tolist():
        elist[t_id][s_id] = year
        rlist[s_id][t_id] = year

    print('gene matrix: ',encoded.shape)
    print('cell matrix: ',encoded2.shape)
    graph.node_feature['gene'] = torch.tensor(encoded, dtype=torch.float)
    graph.node_feature['cell'] = torch.tensor(encoded2, dtype=torch.float)

    graph.years = np.ones(gene_cell.shape[0]+gene_cell.shape[1])
    return graph

def build_data(adj, encoded, encoded2):
    node_type = [0]*adj.shape[0]+[1]*adj.shape[1]
    node_type = torch.LongTensor(node_type)

    g_index,c_index = np.nonzero(adj)
    c_index += adj.shape[0]
    edge_index = torch.tensor([g_index, c_index], dtype=torch.long)
    edge_type = torch.LongTensor([0]*edge_index.shape[1])
    edge_time = torch.LongTensor([0]*edge_index.shape[1])
    
    x = {'gene': torch.tensor(encoded, dtype=torch.float),
         'cell': torch.tensor(encoded2, dtype=torch.float)} 
    # print(len(x['gene']))  # 5000
    # print(len(x['cell']))  # 2713
    # print(len(node_type))  # 7713

    return x,node_type, edge_time, edge_index,edge_type


