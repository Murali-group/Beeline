#!/bin/bash

# Downloading files for lymph node lymnphoma 14k and hg38 and jaspar database into multiomics inputs

if [ ! -e lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5 ]
then
	wget https://bmbl.bmi.osumc.edu/downloadFiles/deepmaps/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5
else
    echo "skipping lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5"
fi

if [ ! -e lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv ]
then
	wget https://bmbl.bmi.osumc.edu/downloadFiles/deepmaps/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv.gz
    gunzip lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv.gz
else
    echo "skipping lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv.gz"
fi

if [ ! -e hg38_1000_2.0.h5 ]
then
	wget http://cistrome.org/~alynch/data/lisa_data/hg38_1000_2.0.h5
else
    echo "skipping hg38_1000_2.h5"
fi

if [ ! -e jaspar_hg38_500.qsave ]
then
	wget https://bmblx.bmi.osumc.edu/downloadFiles/deepmaps/jaspar_hg38_500.qsave
else
    echo "skipping jaspar_hg38_500"
fi

# LISA_PATH=[your_lisa_path]
# $ cd $LISA_PATH
