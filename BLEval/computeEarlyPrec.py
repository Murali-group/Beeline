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

def EarlyPrec(evalObject, algo_name):
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
        trueEdges = set()
        for key in TrueEdgeDict.keys():
            if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                   (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
                    TrueEdgeDict[key] = 1
                    trueEdges.add(key)
                    numEdges += 1

        outDir = str(evalObject.output_settings.base_dir) + \
                 str(evalObject.input_settings.datadir).split("inputs")[1] + \
                 "/" + dataset["name"] + "/" + algo_name

        #algos = evalObject.input_settings.algorithms
        rank_path = outDir + "/rankedEdges.csv"
        if not os.path.isdir(outDir):
            rankDict[dataset["name"]] = set([])
            continue
        try:
            predDF = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
        except:
            print("\nSkipping early precision computation for ", algo_name, "on path", outDir)
            rankDict[dataset["name"]] = set([])
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

    return(np.median(list(Eprec.values())),Eprec)
