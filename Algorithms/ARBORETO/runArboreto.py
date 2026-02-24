from optparse import OptionParser
import os
import sys
import pandas as pd
from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names
from distributed import Client, LocalCluster

def parseArgs(args):
    parser = OptionParser()

    parser.add_option('', '--algo', type = 'str',
                      help='Algorithm to run. Can either by GENIE3 or GRNBoost2')

    parser.add_option('', '--inFile', type='str',
                      help='Path to input tab-separated expression SamplesxGenes file')

    parser.add_option('', '--outFile', type = 'str',
                      help='File where the output network is stored')

    (opts, args) = parser.parse_args(args)

    return opts, args

def main(args):
    opts, args = parseArgs(args)
    inDF = pd.read_csv(opts.inFile, sep = '\t', index_col = 0, header = 0)

    client = Client(processes = False)    

    if opts.algo == 'GENIE3':
        network = genie3(inDF.to_numpy(), client_or_address = client, gene_names = inDF.columns)
        network.to_csv(opts.outFile, index = False, sep = '\t')

    elif opts.algo == 'GRNBoost2':
        network = grnboost2(inDF.to_numpy(), client_or_address = client, gene_names = inDF.columns)
        network.to_csv(opts.outFile, index = False, sep = '\t')

    else:
        print("Wrong algorithm name. Should either be GENIE3 or GRNBoost2.")
                        
if __name__ == "__main__":
    main(sys.argv)
