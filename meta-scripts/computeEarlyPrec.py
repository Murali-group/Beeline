import os
import argparse
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import multiprocessing
from pathlib import Path
import concurrent.futures
from collections import defaultdict
from itertools import product, permutations
from multiprocessing import Pool, cpu_count
from networkx.convert_matrix import from_pandas_adjacency

def EarlyPrec(evalObject, algorithmName, TFEdges = False):
    '''
    Computes early precision for a given algorithm for each dataset.
    We define early precision as the fraction of true 
    positives in the top-k edges, where k is the number of
    edges in the ground truth network (excluding self loops).
    
    
    :param evalObject: An object of class :class:`BLEval.BLEval`.
    :type evalObject: BLEval
      
    :param algorithmName: Name of the algorithm for which the early precision is computed.
    :type algorithmName: str
      
            
    :returns:
        A dataframe containing early precision values
        for a given algorithm for each dataset.

    '''
    rankDict = {}
    for dataset in tqdm(evalObject.input_settings.datasets):
        trueEdgesDF = pd.read_csv(str(evalObject.input_settings.datadir)+'/'+ \
                      dataset['name'] + '/' +\
                      dataset['trueEdges'], sep = ',',
                      header = 0, index_col = None)
        trueEdgesDF = trueEdgesDF.loc[(trueEdgesDF['Gene1'] != trueEdgesDF['Gene2'])]
        trueEdgesDF.drop_duplicates(keep = 'first', inplace=True)
        trueEdgesDF.reset_index(drop=True, inplace=True)


        outDir = str(evalObject.output_settings.base_dir) + \
                 str(evalObject.input_settings.datadir).split("inputs")[1] + \
                 "/" + dataset["name"] + "/" + algorithmName

        #algos = evalObject.input_settings.algorithms
        rank_path = outDir + "/rankedEdges.csv"
        if not os.path.isdir(outDir):
            rankDict[dataset["name"]] = set([])
            continue
        try:
            predDF = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
        except:
            print("\nSkipping early precision computation for ", algorithmName, "on path", outDir)
            rankDict[dataset["name"]] = set([])
            continue

        predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
        predDF.drop_duplicates(keep = 'first', inplace=True)
        predDF.reset_index(drop=True, inplace=True)
        
        
        if TFEdges:
            # Consider only edges going out of TFs
            
            # Get a list of all possible TF to gene interactions 
            uniqueNodes = np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']])
            possibleEdges_TF = set(product(set(trueEdgesDF.Gene1),set(uniqueNodes)))

            # Get a list of all possible interactions 
            possibleEdges_noSelf = set(permutations(uniqueNodes, r = 2))
            
            # Find intersection of above lists to ignore self edges
            # TODO: is there a better way of doing this?
            possibleEdges = possibleEdges_TF.intersection(possibleEdges_noSelf)
            
            TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}

            trueEdges = trueEdgesDF['Gene1'] + "|" + trueEdgesDF['Gene2']
            trueEdges = trueEdges[trueEdges.isin(TrueEdgeDict)]
            print("\nEdges considered ", len(trueEdges))
            numEdges = len(trueEdges)
        
            predDF['Edges'] = predDF['Gene1'] + "|" + predDF['Gene2']
            # limit the predicted edges to the genes that are in the ground truth
            predDF = predDF[predDF['Edges'].isin(TrueEdgeDict)]

        else:
            trueEdges = trueEdgesDF['Gene1'] + "|" + trueEdgesDF['Gene2']
            trueEdges = set(trueEdges.values)
            numEdges = len(trueEdges)
        
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
            print("\nSkipping early precision computation for on path ", rank_path,"due to lack of predictions.")
            rankDict[dataset["name"]] = set([])
    Eprec = {}
    Erec = {}
    for dataset in tqdm(evalObject.input_settings.datasets):
        if len(rankDict[dataset["name"]]) != 0:
            intersectionSet = rankDict[dataset["name"]].intersection(trueEdges)
            Eprec[dataset["name"]] = len(intersectionSet)/len(rankDict[dataset["name"]])
            Erec[dataset["name"]] = len(intersectionSet)/len(trueEdges)
        else:
            Eprec[dataset["name"]] = 0
            Erec[dataset["name"]] = 0

    return(Eprec)
