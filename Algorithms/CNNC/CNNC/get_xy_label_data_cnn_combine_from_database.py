# Usage: python get_xy_data_cnn_combine_from_database.py bulk_gene_list.txt sc_gene_list.txt gene_pair_list  data_separation_index_list  bulk_expression_data  sc_exprsssion_data flag (0, no label. 1, label)
# command line in developer's linux machine :
# python get_xy_label_data_cnn_combine_from_database.py bulk_gene_list.txt sc_gene_list.txt mmukegg_new_new_unique_rand_labelx_sy.txt mmukegg_new_new_unique_rand_labelx_num_sy.txt /home/yey3/sc_process_1/new_bulk_mouse/prs_calculation/mouse_bulk.h5 /home/yey3/sc_process_1/rank_total_gene_rpkm.h5 1
#################INPUT################################################################################################################################
# 1, bulk_gene_list.txt is the list that convert bulk expression data gene set into gene symbol IDs. format: 'gene symbol IDs\t bulk gene ID'
# 2, sc_gene_list.txt is the list that convert sc expression data gene set into gene symbol IDs. format: 'gene symbol IDs\t sc gene ID'
# 3, gene_pair_list is the list that contains gene pairs and their labels. format : 'GeneA    GeneB     0'
# 4, data_separation index list is a number list that divide gene_pair_list into small parts
# here we use data separation index list to divide gene pairs into small data parts, and make sure that the gene pairs in each index inteval is completely isolated from others. And we can evaluate CNNC's performance on only a small data part.
# if users do not need to separate data, they can just generate a index list to divide the data into N equal parts.
# 5, bulk expression data  it should be a hdf5 format. users can use their own data or data we provided.
# 6, sc expression data  it should be a hdf5 format. users can use their own data or data we provided.
# 7, flag (0, no label. 1, label)
#################OUTPUT
# it generate a data_label folder, and a series of data files containing Nxdata_tf (NEPDF file), ydata_tf (label file) and zdata_tf (gene symbol pair file) for each data part divided.

import pandas as pd
from numpy import *
import json, re,os, sys
save_dir = os.path.join(os.getcwd(),'NEPDF_data')
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)
def get_gene_list_bulk(file_name):
    import re
    h={}
    s = open(file_name,'r')   #gene symbol ID list of bulk RNA-seq
    for line in s:
        search_result = re.search(r'^([^\s]+)\s+([^\s]+)',line)
        h[search_result.group(1).lower()]=search_result.group(2)   # h [gene symbol] = gene ID
    s.close()
    return h

def get_gene_list(file_name):
    import re
    h={}
    s = open(file_name,'r') #gene symbol ID list of sc RNA-seq
    for line in s:
        search_result = re.search(r'^([^\s]+)\s+([^\s]+)',line)
        h[search_result.group(1).lower()]=search_result.group(2) # h [gene symbol] = gene ID
    s.close()
    return h

def get_sepration_index (file_name):
    import numpy as np
    index_list = []
    s = open(file_name, 'r')
    for line in s:
        index_list.append(int(line))
    return (np.array(index_list))

# Script starts from here
if len(sys.argv) < 8:
    print ('No enough input files')
    sys.exit()
if not sys.argv[1] == 'None':
    h_gene_list_bulk =get_gene_list_bulk(sys.argv[1]) #'bulk_gene_list.txt')#
    print ('read bulk gene list')
elif sys.argv[5] == 'None':  ### bulk list = none
    print ('No bulk gene list')
else :
    print('wrong bulk expression information')
    sys.exit()
if not sys.argv[2] == 'None':
    h_gene_list =get_gene_list(sys.argv[2]) # 'sc_gene_list.txt')#
    print ('read sc gene list')
elif sys.argv[6] == 'None':  ### sc list = none
    print('No sc gene list')
else:
    print('wrong sc expression information')
    sys.exit()

if not sys.argv[5] == 'None':
    store = pd.HDFStore(sys.argv[5])#'/home/yey3/sc_process_1/new_bulk_mouse/prs_calculation/mouse_bulk.h5')  ### bulk RNA-seq expression data        )#
    rpkm_bulk = store['rpkm']
    store.close()
    print('read bulk RNA-seq expression')
elif sys.argv[1] == 'None':  ### sc list = none
    print('No bulk gene expression')
else:
    print('wrong bulk expression information')
    sys.exit()
if not sys.argv[6] == 'None':
    store = pd.HDFStore(sys.argv[6])#'/home/yey3/sc_process_1/rank_total_gene_rpkm.h5')    # scRNA-seq expression data                        )#
    rpkm = store['rpkm']
    store.close()
    print('read sc RNA-seq expression')
elif sys.argv[2] == 'None':  ### sc list = none
    print('No sc gene expression ')
else:
    print('wrong sc expression information')
    sys.exit()

if sys.argv[1] == 'None' and sys.argv[2] == 'None':
    print ('no bulk or sc data')
    sys.exit()
if sys.argv[5] == 'None' and sys.argv[6] == 'None':
    print ('no bulk or sc data')
    sys.exit()
########## generate NEPDF matrix
gene_pair_label = []
s=open(sys.argv[3])#'mmukegg_new_new_unique_rand_labelx.txt')#)   ### read the gene pair and label file
for line in s:
    gene_pair_label.append(line)
gene_pair_index = get_sepration_index(sys.argv[4])#'mmukegg_new_new_unique_rand_labelx_num.npy')#sys.argv[6]) # read file speration index
s.close()
gene_pair_label_array = array(gene_pair_label)
for i in range(len(gene_pair_index)-1):   #### many sperations
    print (i)
    start_index = gene_pair_index[i]
    end_index = gene_pair_index[i+1]
    x = []
    y = []
    z = []
    for gene_pair in gene_pair_label_array[start_index:end_index]: ## each speration
        separation = gene_pair.split()
        if sys.argv[7] == '1':
            x_gene_name,y_gene_name,label = separation[0],separation[1],separation[2]
            y.append(label)
        else:
            x_gene_name, y_gene_name = separation[0], separation[1]
        z.append(x_gene_name+'\t'+y_gene_name)

        if not sys.argv[1] == 'None':
            x_tf_bulk = log10(rpkm_bulk[h_gene_list_bulk[x_gene_name]][0:249] + 10 ** -2)  ## 249 means the number of samples, users can just remove '[0:249]'
            x_gene_bulk = log10(rpkm_bulk[h_gene_list_bulk[y_gene_name]][0:249] + 10 ** -2)
            H_T_bulk = histogram2d(x_tf_bulk, x_gene_bulk, bins=32)
            H_bulk= H_T_bulk[0].T
            HT_bulk = (log10(H_bulk / 43261 + 10 ** -4) + 4)/4
        if not sys.argv[2] == 'None':
            x_tf = log10(rpkm[int(h_gene_list[x_gene_name])][0:43261] + 10 ** -2) # ## 43261 means the number of samples in the sc data, we also have one row that is sum of all cells, so the real size is 43262, that is why we use [0:43261]. For TF target prediction or other data, just remove "[0:43261]"
            x_gene = log10(rpkm[int(h_gene_list[y_gene_name])][0:43261] + 10 ** -2)# For TF target prediction, remove "[0:43261]"
            H_T = histogram2d(x_tf, x_gene, bins=32)
            H = H_T[0].T
            HT = (log10(H / 43261 + 10 ** -4) + 4) / 4
        if sys.argv[1] == 'None': ## bulk is none, only sc
            x.append(HT)
        elif sys.argv[2] == 'None':  ## sc is none, only bulk
            x.append(HT_bulk)
        else:
            x.append(concatenate((HT, HT_bulk), axis=0))
    if (len(x)>0):
        xx = array(x)[:, :, :, newaxis]
    else:
        xx = array(x)
    save(save_dir+'/Nxdata_tf' + str(i) + '.npy', xx)
    if sys.argv[7] == '1':
        save(save_dir+'/ydata_tf' + str(i) + '.npy', array(y))
    save(save_dir+'/zdata_tf' + str(i) + '.npy', array(z))



