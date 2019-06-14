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
import src.evalClass as ev

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

    parser.add_argument('-t', '--time', action="store_true", default=False,
      help="Analyze time taken by each algorithm for a.\n")
    
    parser.add_argument('-e', '--EPR', action="store_true", default=False,
      help="Compute median early precision.")
    
    parser.add_argument('-s','--sEPR', action="store_true", default=False,
      help="Analyze median (signed) early precision for activation and inhibitory edges.")

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
        
    print('Post-run evaluation started...\n')
    evalSummarizer = ev.Evaluation(evalConfig.input_settings, evalConfig.output_settings)
    
    outDir = str(evalSummarizer.output_settings.base_dir) + \
            str(evalSummarizer.input_settings.datadir).split("inputs")[1] + "/"+\
            str(evalSummarizer.output_settings.output_prefix) + "-"
    
    # Compute and plot ROC, PRC and report median AUROC, AUPRC    
    if (opts.auc):
        print('Computing areas under ROC and PR curves...\n')

        AUPRCDict, AUROCDict, uAUPRCDict, uAUROCDict = evalSummarizer.findAUC()
        pd.DataFrame(AUPRCDict).to_csv(outDir+'AUPRCscores.csv')
        pd.DataFrame(AUROCDict).to_csv(outDir+'AUROCscores.csv')
    
    # Compute Jaccard index    
    if (opts.jaccard):
        print('Computing Jaccard index...\n')

        jaccDict = evalSummarizer.findJaccard()
        jaccDict.to_csv(outDir + "Jaccard.csv")
        
    # Compute Spearman correlation scores
    if (opts.spearman):
        print('Computing Spearman\'s correlation...\n')

        corrDict = evalSummarizer.findSpearman()
        corrDict.to_csv(outDir + "Spearman.csv")
        
    # Compute median time taken
    if (opts.time):
        print('Computing time taken...\n')

        TimeDict = evalSummarizer.parseTimes()
        pd.DataFrame(TimeDict).to_csv(outDir+'Times.csv')
    
    # Compute early precision
    if (opts.EPR):
        print('Computing early precision values...\n')
        ePRDF = evalSummarizer.findEarlyPr()
        ePRDF.to_csv(outDir + "EPr.csv")
                        
    # Compute early precision for activation and inhibitory edges
    if (opts.sEPR):
        print('Computing early precision values for activation and inhibitory edges...\n')
        
        sPrDF = evalSummarizer.findSignedPr()
        sPrDF.to_csv(outDir + "sEPr.csv")

    print('Evaluation complete...\n')


if __name__ == '__main__':
  main()
