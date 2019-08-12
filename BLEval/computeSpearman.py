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

def Spearman(evalObject, algorithmName):
    """
    A function to compute median pairwirse Spearman correlation
    of predicted ranked edges, i.e., the outputs of different datasets
    generated from the same reference network, for a given algorithm.
    

    :param evalObject: An object of class :class:`BLEval.BLEval`.
    :type evalObject: BLEval
      
    :param algorithmName: Name of the algorithm for which the Spearman correlation is computed.
    :type algorithmName: str
      
      
    :returns:
        - median: Median of Spearman correlation values
        - mad: Median Absolute Deviation of  the Spearman correlation values
    """

    rankDict = {}
    sim_names = []
    for dataset in tqdm(evalObject.input_settings.datasets):
        trueEdgesDF = pd.read_csv(str(evalObject.input_settings.datadir)+'/'+ \
                      dataset['name'] + '/' +\
                      dataset['trueEdges'], sep = ',',
                      header = 0, index_col = None)
        possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                     r = 2))
        PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}

        outDir = str(evalObject.output_settings.base_dir) + \
                 str(evalObject.input_settings.datadir).split("inputs")[1] + \
                 "/" + dataset["name"] + "/" + algorithmName
        #algos = evalObject.input_settings.algorithms
        rank_path = outDir+"/rankedEdges.csv"
        if not os.path.isdir(outDir):
            continue
        try:
            predEdgeDF = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
        except:
            print("Skipping spearman computation for ", algorithmName, "on path", outDir)
            continue

        for key in PredEdgeDict.keys():
            subDF = predEdgeDF.loc[(predEdgeDF['Gene1'] == key.split('|')[0]) &
                           (predEdgeDF['Gene2'] == key.split('|')[1])]
            if len(subDF)>0:
                PredEdgeDict[key] = np.abs(subDF.EdgeWeight.values[0])
        rankDict[dataset["name"]] = PredEdgeDict
        sim_names.append(dataset["name"])
    df2 = pd.DataFrame.from_dict(rankDict)
    spearmanDF = df2.corr(method='spearman')


    df = spearmanDF.where(np.triu(np.ones(spearmanDF.shape),  k = 1).astype(np.bool))

    df = df.stack().reset_index()
    df.columns = ['Row','Column','Value']

    return(df.Value.median(),df.Value.mad())
