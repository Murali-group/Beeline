from .conv import *

class Classifier(nn.Module):
    def __init__(self, n_hid, n_out):
        super(Classifier, self).__init__()
        self.n_hid    = n_hid
        self.n_out    = n_out
        self.linear   = nn.Linear(n_hid,  n_out)
    def forward(self, x):
        tx = self.linear(x)
        return torch.log_softmax(tx.squeeze(), dim=-1)
    def __repr__(self):
        return '{}(n_hid={}, n_out={})'.format(
            self.__class__.__name__, self.n_hid, self.n_out)

class Matcher(nn.Module):
    '''
        Matching between a pair of nodes to conduct link prediction.
        Use multi-head attention as matching model.
    '''
    def __init__(self, n_hid):
        super(Matcher, self).__init__()
        self.left_linear    = nn.Linear(n_hid,  n_hid)
        self.right_linear   = nn.Linear(n_hid,  n_hid)
        self.sqrt_hd  = math.sqrt(n_hid)
        self.cache      = None
    def forward(self, x, y, infer = False, pair = False):
        ty = self.right_linear(y)
        if infer:
            '''
                During testing, we will consider millions or even billions of nodes as candidates (x).
                It's not possible to calculate them again for different query (y)
                Since the model is fixed, we propose to cache them, and dirrectly use the results.
            '''
            if self.cache != None:
                tx = self.cache
            else:
                tx = self.left_linear(x)
                self.cache = tx
        else:
            tx = self.left_linear(x)
        if pair:
            res = (tx * ty).sum(dim=-1)
        else:
            res = torch.matmul(tx, ty.transpose(0,1))
        return res / self.sqrt_hd
    def __repr__(self):
        return '{}(n_hid={})'.format(
            self.__class__.__name__, self.n_hid)
    


        
class GNN(nn.Module):
    def __init__(self, in_dim, n_hid, num_types, num_relations, n_heads, n_layers, dropout = 0.2, conv_name = 'hgt', prev_norm = True, last_norm = True, use_RTE = True):
        super(GNN, self).__init__()
        self.gcs = nn.ModuleList()
        self.num_types = num_types
        self.in_dim    = in_dim
        self.n_hid     = n_hid
        self.adapt_ws  = nn.ModuleList()
        self.drop      = nn.Dropout(dropout)
        for t in range(num_types):
            self.adapt_ws.append(nn.Linear(in_dim, n_hid))
        for l in range(n_layers - 1):
            self.gcs.append(GeneralConv(conv_name, n_hid, n_hid, num_types, num_relations, n_heads, dropout, use_norm = prev_norm, use_RTE = use_RTE))
        self.gcs.append(GeneralConv(conv_name, n_hid, n_hid, num_types, num_relations, n_heads, dropout, use_norm = last_norm, use_RTE = use_RTE))

    def forward(self, node_feature, node_type, edge_time, edge_index, edge_type):
        res = torch.zeros(node_feature.size(0), self.n_hid).to(node_feature.device)
        for t_id in range(self.num_types):
            idx = (node_type == int(t_id))
            if idx.sum() == 0:
                continue
            res[idx] = torch.tanh(self.adapt_ws[t_id](node_feature[idx]))
        meta_xs = self.drop(res)
        del res
        for gc in self.gcs:
            meta_xs = gc(meta_xs, node_type, edge_index, edge_type, edge_time)
        return meta_xs  

class GNN_from_raw(nn.Module):
    def __init__(self, in_dim, n_hid, num_types, num_relations, n_heads, n_layers, dropout = 0.2, conv_name = 'hgt', prev_norm = True, last_norm = True, use_RTE = True):
        super(GNN_from_raw, self).__init__()
        self.gcs = nn.ModuleList()
        self.num_types = num_types
        self.in_dim    = in_dim
        self.n_hid     = n_hid
        self.adapt_ws  = nn.ModuleList()
        self.drop      = nn.Dropout(dropout)
        self.embedding1 = nn.ModuleList()
        self.embedding2 = nn.ModuleList()
        #self.fc1 = nn.Linear(dim, 512)
        #self.fc2 = nn.Linear(512, 256)
        #print("in_dim in GNN_from_raw\n")
        #print(in_dim[0]) #2713
        #print(in_dim[1]) #24022
        for ti in range(num_types):
             #self.embedding.append(F.relu(nn.Linear(512,256)(F.relu(nn.Linear(in_dim[ti],512)))))
             self.embedding1.append(nn.Linear(in_dim[ti],512)) #embedding1[0] [2713 x 512] embedding1[1] [24022 x 512]
             self.embedding2.append(nn.Linear(512,256)) #embedding2[0] [512, 256] embedding2[1] [512,256]
        
        for t in range(num_types):
            self.adapt_ws.append(nn.Linear(256, n_hid)) #256 could be one additional hyperparameter!!!
        for l in range(n_layers - 1):
            self.gcs.append(GeneralConv(conv_name, n_hid, n_hid, num_types, num_relations, n_heads, dropout, use_norm = prev_norm, use_RTE = use_RTE))
        self.gcs.append(GeneralConv(conv_name, n_hid, n_hid, num_types, num_relations, n_heads, dropout, use_norm = last_norm, use_RTE = use_RTE))
    
    
    def encode(self, x,t_id):
        h1 = F.relu(self.embedding1[t_id](x))
        return F.relu(self.embedding2[t_id](h1))
    
    def forward(self, node_feature, node_type, edge_time, edge_index, edge_type):
        node_embedding=[] #len = 2 
        for t_id in range(self.num_types):
            node_embedding += list(self.encode(node_feature[t_id],t_id))
        
        node_embedding = torch.stack(node_embedding)
        #print("shape of node_embedding="+str(node_embedding.shape)+"\n")
        res = torch.zeros(node_embedding.size(0), self.n_hid).to(node_feature[0].device)
        
        for t_id in range(self.num_types):
            idx = (node_type == int(t_id)) #0, 1
            if idx.sum() == 0:
                continue
            #res[idx] = torch.tanh(self.adapt_ws[t_id](self.embedding[t_id](node_feature[idx])))
            #print(idx)
            res[idx] = torch.tanh(self.adapt_ws[t_id](node_embedding[idx]))
            #res[idx] = torch.tanh(self.adapt_ws[t_id](self.encode(node_feature[t_id],t_id)))
        
        meta_xs = self.drop(res)
        del res
        for gc in self.gcs:
            meta_xs = gc(meta_xs, node_type, edge_index, edge_type, edge_time)
        return meta_xs  
