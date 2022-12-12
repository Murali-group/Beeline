# Please refer to https://github.com/SPLab-aviyente/scSGL/blob/main/notebooks/demo.ipynb

import os
import argparse
import pandas as pd #to load read GSD dataset
import numpy as np
import sys
sys.path.append('scSGL') #to add a path to search for the requested module

from pysrc.graphlearning import learn_signed_graph
from pysrc.evaluation import auc #to evaluate inference with auprc/auroc


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Run scSGL algorithm.') 
    parser.add_argument('--expression_file', 
        help='Path to ExpressionData file')
    parser.add_argument('--ref_net_file', 
        help='Path to refNetwork file')
    parser.add_argument('--pos_density', default='0.45', #to control the density of positive part of the learned signed graph
        help='Positive density')
    parser.add_argument('--neg_density', default='0.45', #to control the density of negative part of the learned signed graph
        help='Negative density')
    parser.add_argument('--assoc', default='correlation', #to infer a signed graph with correlation kernel
        help='Association type')
    parser.add_argument('--out_file',
        help='Path to output file')
    
    return parser

def parse_arguments():
    parser = get_parser()
    opts = parser.parse_args()

    return opts

def main(args):
    #python run_scSGL.py --expression_file scSGL/data/inputs/GSD/ExpressionData.csv --ref_net_file scSGL/data/inputs/GSD/refNetwork.csv --out_file outFile.txt

    opts = parse_arguments()
    expression_df = pd.read_csv(opts.expression_file, index_col=0)  #to read gene expression file
    ref_net_df = pd.read_csv(opts.ref_net_file) #to read reference network file

    #Learn signed graph with the parameters
    G = learn_signed_graph(expression_df.to_numpy(), pos_density=float(opts.pos_density), neg_density=float(opts.neg_density),
                                assoc=opts.assoc, gene_names=np.array(expression_df.index))
    #G is a dataframe with each row indicating an edge between two genes. 
    #Each edge is also associated with a weight, which is either positive or negative depending on the sign of the edge.
    
    G.to_csv(opts.out_file, index = False, sep = '\t')  #to write the output file

if __name__ == "__main__":
    main(sys.argv)
