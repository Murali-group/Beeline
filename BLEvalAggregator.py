#!/usr/bin/env python
# coding: utf-8

import os
import yaml
import argparse
import itertools
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import multiprocessing
from pathlib import Path
import concurrent.futures
from itertools import permutations
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from networkx.convert_matrix import from_pandas_adjacency

# local imports
import BLEval as ev

def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('-c','--config', default='config.yaml',
        help="Configuration file containing list of datasets "
              "algorithms and output specifications.\n")

    parser.add_argument('-a', '--auc', action="store_true", default=False,
        help="Compute median of areas under Precision-Recall and ROC curves.\n")

    parser.add_argument('-j', '--jaccard', action="store_true", default=False,
      help="Compute median Jaccard index of predicted top-k networks "
      "for each algorithm for a given set of datasets generated "
      "from the same ground truth network.\n")

    parser.add_argument('-r', '--spearman', action="store_true", default=False,
      help="Compute median Spearman Corr. of predicted edges "
      "for each algorithm  for a given set of datasets generated "
      " from the same ground truth network.\n")

    parser.add_argument('-l', '--clr', action="store_true", default=False,
      help="Compute the infered ranked list using the CLR network algorithm"
      "for each algorithm for a given set of datasets generated "
      " from the same ground truth network.\n")

    parser.add_argument('-t', '--time', action="store_true", default=False,
      help="Analyze time taken by each algorithm for a.\n")

    parser.add_argument('-e', '--EPR', action="store_true", default=False,
      help="Compute median early precision.")

    parser.add_argument('-s','--sEPR', action="store_true", default=False,
      help="Analyze median (signed) early precision for activation and inhibitory edges.")

    parser.add_argument('-m','--motifs', action="store_true", default=False,
      help="Compute network motifs in the predicted top-k networks.")

    parser.add_argument('-b','--borda', action="store_true", default=False,
      help="Compute edge ranked list using the borda method.")

    parser.add_argument('--borda-agg', action="store", default="average",
      help="Method used to aggregate rank in borda method. Available options are {‘average’, ‘min’, ‘max’, ‘first’}, default ‘average’")

    parser.add_argument('--borda-algo', action="append",
      help="Algorithms used to compute rank in borda method. Default option runs Borda on all algorithm.")

    parser.add_argument('-i', '--indirect-auc', action="store_true", default=False,
        help="Compute areas under k-Precision-Recall and k-ROC curves while accounting for indirect paths of length k as true edges.\n")

    parser.add_argument('-k', '--indirect-auc-k', action="store", default=1,
        help="A predicted edge (a,b) is considered a TP if there is a path of length less than equal to k between a and b. So a 1-PR curve is just the regular PR curve.\n")

    parser.add_argument('-R','--ref', default=None,
        help="Reference network file containing list of true positive (TP) edges. This option can be used in conjuction with `--auc` or `--indirect-auc` option to compute AUPRC and AUROC scores with given network as reference. The file should be comma separated and have following node column names: `Gene1` and `Gene2`.\n")

    parser.add_argument('-u','--undirected', action="store_true", default=False,
      help="A flag to indicate whether to treat predictions as undirected edges (undirected = True) or directed edges (undirected = False) during AUC computations. This option can be used in conjuction with `--auc` or `--indirect-auc`. ")

    return parser

def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts

def main():
    opts = parse_arguments()
    config_file = opts.config

    evalConfig = None

    with open(config_file, 'r') as conf:
        evalConfig = ev.ConfigParser.parse(conf)

    print('\nPost-run evaluation started...')
    evalSummarizer = ev.BLEval(evalConfig.input_settings, evalConfig.output_settings)

    outDir = str(evalSummarizer.output_settings.base_dir) + \
            str(evalSummarizer.input_settings.datadir).split("inputs")[1] + "/"+\
            str(evalSummarizer.output_settings.output_prefix) + "-"

    # Compute and plot ROC, PRC and report median AUROC, AUPRC
    if (opts.auc):
        print('\n\nComputing areas under ROC and PR curves...')

        AUPRC, AUROC = evalSummarizer.computeAUC(userReferenceNetworkFile=opts.ref, directed=not opts.undirected)
        AUPRC.to_csv(outDir+'AUPRC.csv')
        AUROC.to_csv(outDir+'AUROC.csv')

    # Compute Jaccard index
    if (opts.jaccard):
        print('\n\nComputing Jaccard index...')

        jaccDict = evalSummarizer.computeJaccard()
        jaccDict.to_csv(outDir + "Jaccard.csv")

    # Compute Spearman correlation scores
    if (opts.spearman):
        print('\n\nComputing Spearman\'s correlation...')

        corrDict = evalSummarizer.computeSpearman()
        corrDict.to_csv(outDir + "Spearman.csv")

    # Compute inferred ranked edge list using CLR algorithm
    if (opts.clr):
        print('\n\nComputing inferred ranked edge list using CLR algorithm...')

        AUPRC, AUROC = evalSummarizer.computeCLR()
        AUPRC.to_csv(outDir + 'CLR-AUPRC.csv')
        AUROC.to_csv(outDir + 'CLR-AUROC.csv')

    # Compute median time taken
    if (opts.time):
        print('\n\nComputing time taken...')

        TimeDict = evalSummarizer.parseTime()
        pd.DataFrame(TimeDict).to_csv(outDir+'Times.csv')

    # Compute early precision
    if (opts.EPR):
        print('\n\nComputing early precision values...')
        ePRDF = evalSummarizer.computeEarlyPrec()
        ePRDF.to_csv(outDir + "EPr.csv")

    # Compute early precision for activation and inhibitory edges
    if (opts.sEPR):
        print('\n\nComputing early precision values for activation and inhibitory edges...')

        actDF, inhDF = evalSummarizer.computeSignedEPrec()
        actDF.to_csv(outDir + "EPr-Activation.csv")
        inhDF.to_csv(outDir + "EPr-Inhibitory.csv")

    # Compute median time taken
    if (opts.motifs):
        print('\n\nComputing network motifs...')

        FBL, FFL, MI = evalSummarizer.computeNetMotifs()
        FBL.to_csv(outDir+'NetworkMotifs-FBL.csv')
        FFL.to_csv(outDir+'NetworkMotifs-FFL.csv')
        MI.to_csv(outDir+'NetworkMotifs-MI.csv')

    # Compute edge ranked list using the borda method
    if (opts.borda):
        print('\n\nComputing edge ranked list using the borda method')
        evalSummarizer.computeBorda(selectedAlgorithms=opts.borda_algo, aggregationMethod=opts.borda_agg)

    # Computing areas under k-Precision-Recall and k-ROC curves while accounting for indirect paths of length k as TPs
    if (opts.indirect_auc):
        print('\n\nComputing areas under %s-Precision-Recall and %s-ROC curves' % (opts.indirect_auc_k, opts.indirect_auc_k))
        AUPRC, AUROC = evalSummarizer.computekAUC(k=int(opts.indirect_auc_k), directed=not opts.undirected, userReferenceNetworkFile=opts.ref)
        AUPRC.to_csv(outDir + '%s-AUPRC.csv' % opts.indirect_auc_k)
        AUROC.to_csv(outDir + '%s-AUROC.csv' % opts.indirect_auc_k)


    print('\n\nEvaluation complete...\n')


if __name__ == '__main__':
  main()
