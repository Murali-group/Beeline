import os
import argparse
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

def Jaccard(evalObject, algorithmName):
    """
    A function to compute median pairwirse Jaccard similarity index
    of predicted top-k edges for a given set of datasets (obtained from
    the same reference network). Here k is the number of edges in the
    reference network (excluding self loops). 
    
    
    :param evalObject: An object of class :class:`BLEval.BLEval`.
    :type evalObject: :obj:`BLEval`
      
      
    :param algorithmName: Name of the algorithm for which the Spearman correlation is computed.
    :type algorithmName: str
      
      
    :returns:
        - median: Median of Jaccard correlation values
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

        TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}

        # Compute TrueEdgeDict Dictionary
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth
        numEdges = 0
        for key in TrueEdgeDict.keys():
            if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                   (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
                    TrueEdgeDict[key] = 1
                    numEdges += 1

        outDir = str(evalObject.output_settings.base_dir) + \
                 str(evalObject.input_settings.datadir).split("inputs")[1] + \
                 "/" + dataset["name"] + "/" + algorithmName

        #algos = evalObject.input_settings.algorithms
        rank_path = outDir + "/rankedEdges.csv"
        if not os.path.isdir(outDir):
            continue
        try:
            predDF = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
        except:
            print("Skipping Jaccard computation for ", algorithmName, "on path", outDir)
            continue

        predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
        predDF.drop_duplicates(keep = 'first', inplace=True)
        predDF.reset_index(drop = True,  inplace= True)
        # check if ranked edges list is empty
        # if so, it is just set to an empty set

        if not predDF.shape[0] == 0:

            # we want to ensure that we do not include
            # edges without any edge weight
            # so check if the non-zero minimum is
            # greater than the edge weight of the top-kth
            # node, else use the non-zero minimum value.
            predDF.EdgeWeight = predDF.EdgeWeight.round(6)
            predDF.EdgeWeight = predDF.EdgeWeight.abs()

            # Use num True edges or the number of
            # edges in the dataframe, which ever is lower
            maxk = min(predDF.shape[0], numEdges)
            edgeWeightTopk = predDF.iloc[maxk-1].EdgeWeight

            nonZeroMin = np.nanmin(predDF.EdgeWeight.replace(0, np.nan).values)
            bestVal = max(nonZeroMin, edgeWeightTopk)

            newDF = predDF.loc[(predDF['EdgeWeight'] >= bestVal)]
            rankDict[dataset["name"]] = set(newDF['Gene1'] + "|" + newDF['Gene2'])
        else:
            rankDict[dataset["name"]] = set([])

    Jdf = computePairwiseJacc(rankDict)
    df = Jdf.where(np.triu(np.ones(Jdf.shape),  k = 1).astype(np.bool))
    df = df.stack().reset_index()
    df.columns = ['Row','Column','Value']
    return(df.Value.median(),df.Value.mad())


def computePairwiseJacc(inDict):
    """
    A helper function to compute all pairwise Jaccard similarity indices
    of predicted top-k edges for a given set of datasets (obtained from
    the same reference network). Here k is the number of edges in the
    reference network (excluding self loops). 
    
    

    :param inDict:  A dictionary contaninig top-k predicted edges  for each dataset. Here, keys are the dataset name and the values are the set of top-k edges.
    :type inDict: dict
    :returns:
        A dataframe containing pairwise Jaccard similarity index values
    """
    jaccDF = {key:{key1:{} for key1 in inDict.keys()} for key in inDict.keys()}
    for key_i in inDict.keys():
        for key_j in inDict.keys():
            num = len(inDict[key_i].intersection(inDict[key_j]))
            den = len(inDict[key_i].union(inDict[key_j]))
            if den != 0:
                jaccDF[key_i][key_j] = num/den
            else:
                jaccDF[key_i][key_j] = 0
    return pd.DataFrame(jaccDF)
