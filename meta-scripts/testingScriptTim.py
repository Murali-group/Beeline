import pandas as pd
import numpy as np
import split_script as split
import os

ExpressionData = pd.read_csv("inputs/human-scRNA/hESC/hESC-500-ChIP/hESC-500-ChIP-ExpressionData.csv",
                                     header = 0, index_col = 0)
gene_pairs = pd.read_csv("inputs/human-scRNA/hESC/hESC-500-ChIP/hESC-500-ChIP-network.csv", header = 0) # check for directed column
# GenePairs = GenePairs.drop(columns='Type', inplace=True)

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
print(final_pairs)
final_pairs.to_csv("./genes.txt", header=False, index=False, sep='\t')

train, test = split.split_edge_cv(final_pairs, 3, 1)
for i in range(3):
    split.verify_split(final_pairs, train[i], test[i], "edge")

# os.system("python 'Algorithms/CNNC/CNNC/get_xy_label_data_cnn_combine_from_database.py' None ")






# GenePairs.to_csv("./genes.csv", header=False, index=False) 

# dendritic = pd.read_csv("Algorithms/CNNC/CNNC/data/dendritic_gene_pairs_200.txt", header=None, delimiter='\t')
# dendritic.columns = ['source', 'target', 'label']
# print(dendritic.groupby('source')['target'].agg([('count', 'size')]))
# print(dendritic)
