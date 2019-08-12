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

def signedEPrec(evalObject, algorithmName):
    '''
    Computes median signed early precision for a given algorithm across all datasets, 
    i.e., the function computes early precision of activation edges and
    early precision for inhibitory edges in the reference network.
    We define early precision of activation edges as the fraction of true 
    positives in the top-ka edges, where ka is the number of activation
    edges in the reference network (excluding self loops). 
    We define early precision of inhibitory edges as the fraction of true 
    positives in the top-ki edges, where ki is the number of inhibitory
    edges in the reference network (excluding self loops).
    

    :param evalObject: An object of class :class:`BLEval.BLEval`.
    :type evalObject: BLEval
      
    :param algorithmName: Name of the algorithm for which the early precision is computed.
    :type algorithmName: str
      
            
    :returns: 
        A dataframe with early precision of activation edges (+) and inhibitory edges (-)
        for a given algorithm
    '''
        
    rankDict = {'+':{},'-':{}}
    sim_names = []
    for sgn in ['+','-']:
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
            trueEdges = set()
            toRemove = []
            for key in TrueEdgeDict.keys():
                subDF = trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                       (trueEdgesDF['Gene2'] == key.split('|')[1])]

                if subDF.shape[0] > 0:
                    if  subDF['Type'].values[0] == sgn:
                        TrueEdgeDict[key] = 1
                        trueEdges.add(key)
                        numEdges += 1
                    else:
                        toRemove.append(key)
            for key in toRemove:
                TrueEdgeDict.pop(key, None)


            outDir = str(evalObject.output_settings.base_dir) + \
                     str(evalObject.input_settings.datadir).split("inputs")[1] + \
                     "/" + dataset["name"] + "/" + algorithmName

            #algos = evalObject.input_settings.algorithms
            rank_path = outDir + "/rankedEdges.csv"
            if not os.path.isdir(outDir):
                print(outDir," not found")
                rankDict[sgn][dataset["name"]] = set([])
                continue
            try:
                predDF = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
            except:
                print("Skipping signed precision computation for ", algorithmName, "on path", outDir)
                rankDict[sgn][dataset["name"]] = set([])
                continue

            predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
            predDF.drop_duplicates(keep = 'first', inplace=True)
            predDF.reset_index(drop = True,  inplace= True)

            # Remove incorrect sign from consideration
            for idx, row in predDF.iterrows():
                if str(row['Gene1']) + '|' + str(row['Gene2']) not in TrueEdgeDict.keys():
                    predDF.drop(idx, axis = 'index', inplace= True)
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
                rankDict[sgn][dataset["name"]] = set(newDF['Gene1'] + "|" + newDF['Gene2'])
            else:
                print("\nSkipping signed early precision computation for file on path ", rank_path,"due to lack of predictions.")
                rankDict[sgn][dataset["name"]] = set([])

    Pprec = {'+':{},'-':{}}
    for sgn in ['+','-']:

        TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}

        # Compute TrueEdgeDict Dictionary
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth
        trueEdges = set()
        for key in TrueEdgeDict.keys():
            subDF = trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                   (trueEdgesDF['Gene2'] == key.split('|')[1])]

            if subDF.shape[0] > 0:
                if  subDF['Type'].values[0] == sgn:
                    TrueEdgeDict[key] = 1
                    trueEdges.add(key)

        for dataset in tqdm(evalObject.input_settings.datasets):
            if len(rankDict[sgn][dataset["name"]]) != 0 and len(trueEdges) != 0:
                intersectionSet = rankDict[sgn][dataset["name"]].intersection(trueEdges)
                Pprec[sgn][dataset["name"]] = len(intersectionSet)/len(rankDict[sgn][dataset["name"]])
            else:
                Pprec[sgn][dataset["name"]] = 0
                
    # To return just the median values, uncomment the line below
    #return(pd.DataFrame(Pprec).median(axis='index').\
    #values[0],pd.DataFrame(Pprec).median(axis='index').values[1])

    return(pd.DataFrame(Pprec))


