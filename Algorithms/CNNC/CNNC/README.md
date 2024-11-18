# CNNC
# Title: Deep learning for inferring gene relationships from single-cell expression data
Originally, CNNC is short for Convolutional neural network co-expression analysis. Co-expression is just one of the tasks CNNC can      do,   but the name (CNNC) derived from it looks very nice.

## see `https://www.pnas.org/content/early/2019/12/09/1911536116` for details.

## date: 2019-07-07


># 1, CNNC
![](https://raw.githubusercontent.com/xiaoyeye/CNNC/master/New%20Bitmap%20Image.bmp)

CNNC aims to infer gene-gene relationships using single cell expression data. For
each gene pair, sc RNA-Seq expression levels are transformed into 32×32 normalized
empirical probability function (NEPDF) matrices. The NEPDF serves as an input to a
convolutional neural network (CNN). The intermediate layer of the CNN can be further
concatenated with input vectors representing Dnase-seq and PWM data. The output
layer can either have a single, three or more values, depending on the application. For
example, for causality inference the output layer contains three probability nodes
where p0 represents the probability that genes a and b are not interacting, p1
encodes the case that gene a regulates gene b, and p2 is the probability that gene b
regulates gene a.
># 2, Pipelines
![](https://raw.githubusercontent.com/xiaoyeye/CNNC/master/pipeline.bmp)

(a) Pipeline for TF-target, KEGG and Reactome edge predictions. Users only need to provide gene-pair candidate list. TF-tatget prediction is cell type specific. Here we provide the model for mESC TF prediction. Please use mESC expression data to generate mESC subset  and then do following NEPDF generation and classification of training and test, and use the big scRNA-seq and bulk data to do pathway tasks. (b) Pipeline for a new task with the expression data we collected. Users need to provide gene-pair candidate list to generate NEPDF list and label list to train and test model. (c) Pipeline for a new task with the expression data users collect. Users need to provide gene-pair candidate list, their own expression data to generate NEPDF list, and label list to train and test model. 
># 3, Data sources
>>## 3.1 scRNA-seq : 
    https://s3.amazonaws.com/mousescexpression/rank_total_gene_rpkm.h5

>>## 3.2 bone marrow drived macrophage scRNA-seq : 
    https://mousescexpression.s3.amazonaws.com/bone_marrow_cell.h5

>>## 3.3 bulk RNA-seq : 
    https://s3.us-east-2.amazonaws.com/mousebulkexprssion/mouse_bulk.h5
>>## 3.4 mesc single cell RNA-seq:
    https://mousescexpression.s3.amazonaws.com/mesc_cell.h5
>>## 3.5 dendritic single cell RNA-seq:
    https://mousescexpression.s3.amazonaws.com/dendritic_cell.h5
># 4, Code environment

>>## Users need to install the latest python and all the modules required by the code.  

Author's environment is python 3.6.3 in a Linux server which is now running Centos 6.5
as the underlying OS and Rocks 6.1.1 as the cluster management revision. 

And Author uses theano as the Keras backend in python. 

Author's GPU is GeForce GTX 1080. If the latest theano does not work, please try some older versions.

Although not necessary, we strongly recommend GPU acceleration and conda management for package, dependency and environment to save time. With conda, the total software, package module installation time in Python should be less than one hour.

># 5, Trained model for 
(see folder for details)
>>## 5.1 KEGG Pathway prediction model

>>## 5.2 Reactome Pathway prediction model

>>## 5.3 GTRD bone TF prediction model

># 6, Train model for a new task

>>## Users can define their own tasks by providing new expression data and (or) new gene pair labels.



># 7, Command lines for Trained model

>>## 7.1 step1, users need to provide gene pair list;

`gene_pair_list` is the list that contains gene pairs and their labels. format : `'GeneA GeneB ' or 'GeneA    GeneB     0'`
such as `mmukegg_new_new_unique_rand_labelx_sy.txt` and `mmukegg_new_new_unique_rand_labelx.txt` in data folder.
users also need to provide `data_separation index_list` which is a number list dividing gene_pair_list into small parts.

Here we use `data separation index list` to divide gene pairs into small data parts, and make sure that the gene pairs in each index inteval is completely isolated from others. We can evaluate CNNC's performance on only a small data part.
If users do not want to specified separate data, they can just generate a index list to divide the data into N equal parts.
>>## 7.2 step2, use `get_xy_label_data_cnn_combine_from_database.py` to get gene pair NEPDF list;

>>>### Usage: python get_xy_data_cnn_combine_from_database.py bulk_gene_list.txt sc_gene_list.txt gene_pair_list  data_separation_index_list  bulk_expression_data  sc_exprsssion_data 0

>>>### command line in author's linux machine :

    python get_xy_label_data_cnn_combine_from_database.py /home/yey/CNNC-master/data/bulk_gene_list.txt /home/yey/CNNC-master/data/sc_gene_list.txt /home/yey/CNNC-master/data/mmukegg_new_new_unique_rand_labelx.txt /home/yey/CNNC-master/data/mmukegg_new_new_unique_rand_labelx_num_sy.txt /home/yey/sc_process_1/new_bulk_mouse/prs_calculation/mouse_bulk.h5 /home/yey/sc_process_1/rank_total_gene_rpkm.h5 0 

#################INPUT################################################################################################################################

>>>#1, `bulk_gene_list.txt` is the list that converts bulk expression data gene set into gene symbol IDs. Format: `'gene symbol IDs\t bulk gene ID'`. Set as `None` if you do not have it.

>>>#2, `sc_gene_list.txt` is the list that converts sc expression data gene set into gene symbol IDs. Format: `'gene symbol IDs\t sc gene ID'`. `please notice that mesc single cell data is not from the big 40k data, so we proivde a new sc_gene_list.txt file for it`.
. Set as `None` if you do not have it.

>>>#3, `gene_pair_list` is the list that contains gene pairs and their labels. format : `'GeneA    GeneB'`

>>>#4, `data_separation index_list` is a number list that divides gene_pair_list into small parts

#Here we use data separation index list to divide gene pairs into small data parts, and make sure that the gene pairs in each index inteval is completely isolated from others. And we can evaluate CNNC's performance on only a small data part.

#if users do not need to separate data, they can just generate a index list to divide the data into N equal parts.

>>>#5, `bulk_expression_data`  it should be a hdf5 format. users can use their own data or data we provided. Set as `None` if you do not have it.

>>>#6, `sc expression data`  it should be a hdf5 format. users can use their own data or data we provided. Set as `None` if you do not have it.

>>>#7， `flag`, 0 means do not generate label list; 1 means to generate label list.

>>> If user does not have bulk (single cell) data, just put `None`s for `bulk_gene_list.txt` (`sc_gene_list.txt`) and `bulk_expression_data` (`sc expression data`).

#################OUTPUT

It generates a NEPDF_data folder, and a series of data files containing `Nxdata_tf` (NEPDF file)  and `zdata_tf` (gene symbol pair file) for each data part divided.

Here we use gene symbol information to align bulk, scRNA-seq and gene pair's gene sets. In our own data, scRNA-seq used entrez ID, bulk RNA-seq used ensembl ID, gene pair list used gene symbol ID, thus we used `bulk_gene_list.txt` and `sc_gene_list.txt` to convert all the IDs to gene symbols. Please also do IDs convertion for bulk and scRNA-seq data if users want to use their own expression data.

>>## 7.3 step3, use `predict_no_y.py` to do prediction;

>>>### Usage: python predict_no_y.py  number_of_separation NEPDF_pathway number_of_categories  model_pathway

>>>### command line in author's linux machine :

    python predict_no_y.py  9 /home/yey/CNNC-master/NEPDF_data  3 /home/yey/CNNC-master/trained_models/KEGG_keras_cnn_trained_model_shallow.h5

(In the models folder are  trained models for  KEGG and Reactome database respectively)

># 8, Command lines for Train new model:

>>## 8.1 step1, users need to provide gene pair candidate list (the same to step1 in trained_model except the flag setting);

`gene_pair_list` is the list that contains gene pairs and their labels. format : `'GeneA    GeneB     0'`
such as `mmukegg_new_new_unique_rand_labelx_sy.txt` and `mmukegg_new_new_unique_rand_labelx.txt` in data folder.
users also need to provide data_separation index_list which is a number list dividing gene_pair_list into small parts

Here we use `data separation index list` to divide gene pairs into small data parts, and make sure that the gene pairs in each index inteval is completely isolated from others. And we can evaluate CNNC's performance on only a small data part.

If users do not need to separate data, they can just generate a index list to divide the data into N equal parts.


>>## 8.2 step2, use `get_xy_label_data_cnn_combine_from_database.py` to get gene pair NEPDF list and their labels (the same to step2 in trained_model except flag setting);

>>>### Usage: python get_xy_data_cnn_combine_from_database.py bulk_gene_list.txt sc_gene_list.txt gene_pair_list  data_separation index list  bulk expression data  sc exprsssion data 1 

>>>### command line in author's linux machine :

    python get_xy_label_data_cnn_combine_from_database.py /home/yey/CNNC-master/data/bulk_gene_list.txt /home/yey/CNNC-master/data/sc_gene_list.txt /home/yey/CNNC-master/data/mmukegg_new_new_unique_rand_labelx.txt /home/yey/CNNC-master/data/mmukegg_new_new_unique_rand_labelx_num_sy.txt /home/yey/sc_process_1/new_bulk_mouse/prs_calculation/mouse_bulk.h5 /home/yey/sc_process_1/rank_total_gene_rpkm.h5 1 

#################INPUT################################################################################################################################

>>>#1, `bulk_gene_list.txt` is the list that converts bulk expression data gene set into gene symbol IDs. Format: `'gene symbol IDs\t bulk gene ID'`. Set as `None` if you do not have it.

>>>#2, `sc_gene_list.txt` is the list that converts sc expression data gene set into gene symbol IDs. Format: `'gene symbol IDs\t sc gene ID'`. Set as `None` if you do not have it.

>>>#3, `gene_pair_list` is the list that contains gene pairs and their labels. format : `'GeneA    GeneB     0'`

>>>#4, `data_separation_index_list` is a number list that divide gene_pair_list into small parts

#Here we use data separation index list to divide gene pairs into small data parts, and make sure that the gene pairs in each index inteval is completely isolated from others. And we can evaluate CNNC's performance on only a small data part.

#if users do not want to separate data, they can just generate a index list to divide the data into N equal parts.

>>>#5,` bulk_expression_data`  it should be a hdf5 format. users can use their own data or data we provided. Set as `None` if you do not have it.

>>>#6, `sc_expression_data`  it should be a hdf5 format. users can use their own data or data we provided. Set as `None` if you do not have it.

>>>#7， `flag`, 0 means do not generate label list; 1 means to generate label list.

>>> If user does not have bulk (single cell) data, just put `None`s for `bulk_gene_list.txt` (`sc_gene_list.txt`) and `bulk_expression_data` (`sc expression data`).
#################OUTPUT

It generate a NEPDF_data folder, and a series of data files containing `Nxdata_tf` (NEPDF file), `ydata_tf` (label file) and `zdata_tf` (gene symbol pair file) for each data part divided.


>>## 8.3 step3, use `train_with_labels_three_foldx.py` to train a new model with three-fold cross validation;

>>>### Usage  python train_with_labels_three_foldx.py number_of_data_parts_divided NEPDF_pathway number_of_categories

>>>### command line in author's linux machine :

    module load cuda-8.0 (it is to use GPU)
    srun -p gpu --gres=gpu:1 -c 2 --mem=20Gb python train_with_labels_three_foldx.py 9 /home/yey/CNNC-master/NEPDF_data 3 > results.txt

#######################OUTPUT

It generates three cross_Validation folder whose name begins with `YYYYY`, in which `keras_cnn_trained_model_shallow.h5` is the trained model


># 9 Attentions:
>>## Please read the readme very carefully and make sure that all files have correct paths, since users may have different computer environments.
>>## All command lines are just used for demo. If you want to run the real data, plz replace "mmukegg_new_new_unique_rand_labelx_num_sy.txt" with "mmukegg_new_new_unique_rand_labelx_num.txt", and replace all "9" with 3057, which is the real number of seperations.

>>## When label list is very large, say more than 100,000 gene pairs, we recommend users to feed a series of small  number_of_data_parts_divided to run the NEPDF generation in parallel.
>>## We are exploring new tasks and new datatypes for CNNC, to be continued...


# Enjoy our CNNC!!
