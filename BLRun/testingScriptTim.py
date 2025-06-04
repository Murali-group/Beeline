import sys
import pandas as pd
import numpy as np
import importlib
import os
split = importlib.import_module("split_script")

def get_gene_list(file_name):
    import re
    h={}
    s = open(file_name,'r') #gene symbol ID list of sc RNA-seq
    for line in s:
        search_result = re.search(r'^([^\s]+)\s+([^\s]+)',line)
        h[search_result.group(1).lower()]=search_result.group(2) # h [gene symbol] = gene ID
    s.close()
    return h

save_dir = "./"
expression_data = pd.read_csv("inputs/human-scRNA/hESC/hESC-500-ChIP/hESC-500-ChIP-ExpressionData.csv",
                                     header = 0, index_col = 0).T
# print(ExpressionData.shape[0])
expression_data.index = range(expression_data.shape[0])
expression_data.to_hdf("./expData.h5", key="expData")
# h5_test = pd.read_hdf("mesc_cell.h5")
# print(h5_test)
# print(ExpressionData)
gene_labels = pd.DataFrame(expression_data.columns)
gene_labels.columns = ['c1']
gene_labels['c2'] = gene_labels['c1']
gene_labels.to_csv("./geneLabels.txt", header=None, index=None, sep='\t')

gene_pairs = pd.read_csv("inputs/human-scRNA/hESC/hESC-500-ChIP/hESC-500-ChIP-network.csv", header = 0) # check for directed column
gene_pairs.columns = ['source', 'target']
# print(GenePairs.head())
# print(GenePairs.shape)
negative_pairs = split.generate_negative_samples_CNNC(df=gene_pairs, seed=1)
negative_pairs['label'] = 0
gene_pairs['label'] = 1
# print(negative_pairs)
reversed_pairs = gene_pairs.copy(deep=True)
reversed_pairs[['source', 'target']] = reversed_pairs[['target', 'source']]
reversed_pairs['label'] = 2
# print(reversed_pairs)
final_pairs = pd.concat([gene_pairs, reversed_pairs, negative_pairs], axis=0).sort_index(kind='merge')
# print(final_pairs)
final_pairs.to_csv("./genes.txt", header=False, index=False, sep='\t')
k = 3
train, test = split.split_edge_cv(gene_pairs, k, 1)
# print(final_pairs.loc[test[0]])
# print(expression_data['WDHD1'])
# print(expression_data['ZNF239'])
# print(final_pairs)

# test.to_csv("./testGenes.txt", header=False, index=False, sep='\n')
# print(final_pairs.loc[test[0]])

# store = pd.HDFStore("rank_total_gene_rpkm.h5")#'/home/yey3/sc_process_1/rank_total_gene_rpkm.h5')    # scRNA-seq expression data
# store2 = pd.HDFStore("mesc_cell.h5")    
# print(store.info())
# print(store2.info())
# print(pd.read_hdf(store).head())
# print(expression_data.head())
# print(pd.read_hdf(store2).head())

# rpkm = store['rpkm']
# print(rpkm.head())
# store.close()
# store2.close()
# print(expression_data[int(h_gene_list[x_gene_name])
# new_gene_list = get_gene_list()
# print(get_gene_list("geneLabels.txt"))



# quit() #gg
# TODO gene pair labels = final pairs
# gene_pair_index = get_sepration_index(sys.argv[4])#'mmukegg_new_new_unique_rand_labelx_num.npy')#sys.argv[6]) # read file speration index
# gene_indexes = pd.read_csv(sys.argv[4]) # TODO fix
for i in range(k):   #TODO fix
    split.verify_split(gene_pairs, train[i], test[i], "edge")
    # start_index = gene_pair_index[i]
    # end_index = gene_pair_index[i+1]
    x = []
    y = []
    z = []
    for index, gene_pair in final_pairs.loc[test[i]].iterrows(): ## each speration
        print(gene_pair)
        x_gene_name, y_gene_name, label = gene_pair[0], gene_pair[1], gene_pair[2]
        # print(x_gene_name)
        y.append(label)
        z.append(x_gene_name+'\t'+y_gene_name)
        x_tf = np.log10(expression_data[x_gene_name] + 10 ** -2) # ## 43261 means the number of samples in the sc data, we also have one row that is sum of all cells, so the real size is 43262, that is why we use [0:43261]. For TF target prediction or other data, just remove "[0:43261]"
        x_gene = np.log10(expression_data[y_gene_name] + 10 ** -2)# For TF target prediction, remove "[0:43261]"
        H_T = np.histogram2d(x_tf, x_gene, bins=32)
        H = H_T[0].T
        HT = (np.log10(H / 43261 + 10 ** -4) + 4) / 4
        x.append(HT)
    if (len(x)>0):
        xx = np.array(x)[:, :, :, np.newaxis]
    else:
        xx = np.array(x)
    np.save(save_dir+'/Nxdata_tf' + str(i) + '.npy', xx)
    np.save(save_dir+'/ydata_tf' + str(i) + '.npy', np.array(y))
    np.save(save_dir+'/zdata_tf' + str(i) + '.npy', np.array(z))


# GenePairs.to_csv("./genes.csv", header=False, index=False) 

# dendritic = pd.read_csv("Algorithms/CNNC/CNNC/data/dendritic_gene_pairs_200.txt", header=None, delimiter='\t')
# dendritic.columns = ['source', 'target', 'label']
# print(dendritic.groupby('source')['target'].agg([('count', 'size')]))
# print(dendritic)
