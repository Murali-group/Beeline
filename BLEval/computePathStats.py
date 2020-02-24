import pandas as pd
import sys
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from itertools import product, permutations, combinations, combinations_with_replacement
from tqdm import tqdm
import networkx as nx

def pathAnalysis(dataDict, inputSettings):
    '''
    Computes "directed","feed-forward", 
    "cascade", and "mutual" motifs.
    '''
    
    # Read file for trueEdges
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', 
                                header = 0, index_col = None)
            
    possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                 r = 2))        
    EdgeDict = {'|'.join(p):0 for p in possibleEdges}

    refGraph = nx.DiGraph()

    for key in EdgeDict.keys():
        u = key.split('|')[0]
        v = key.split('|')[1]
        if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == u) &
               (trueEdgesDF['Gene2'] == v)])>0:
                refGraph.add_edge(u,v)

    #refCC, refFB, refFF, refMI = getNetProp(refGraph)

    
    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    #print(dataDict['name'])

    ##################################################
    # Get counts of tp,fp with, and fp without paths #
    ##################################################
    collection = {}

    for algo in inputSettings.algorithms:
        # check if the output rankedEdges file exists
        if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists() and algo[0] not in ['PPCOR','PIDC']:
            # Initialize Precsion
            predDF = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                        sep = '\t', header =  0, index_col = None)

            predDF.EdgeWeight = predDF.EdgeWeight.round(6)
            predDF.EdgeWeight = predDF.EdgeWeight.abs()
            predDF = predDF.loc[(predDF['EdgeWeight'] > 0)]
            predDF.drop_duplicates(keep = 'first', inplace=True)
            predDF.reset_index(drop = True,  inplace= True)
            predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
            if predDF.shape[0] != 0:
                maxk = min(predDF.shape[0], len(refGraph.edges()))
                edgeWeightTopk = predDF.iloc[maxk-1].EdgeWeight            
             
            newDF = predDF.loc[(predDF['EdgeWeight'] >= edgeWeightTopk)]
                
            predGraph = nx.DiGraph()
            

            for key in EdgeDict.keys():
                u = key.split('|')[0]
                v = key.split('|')[1]
                if len(newDF.loc[(newDF['Gene1'] == u) &
                       (newDF['Gene2'] == v)])>0:
                        predGraph.add_edge(u,v)

            dataDict = pathStats(predGraph, refGraph)
            collection[algo[0]] = dataDict
            if algo[0] == 'PPCOR' or algo[0] == 'PIDC':
                collection[algo[0]] = {}
        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' is an undirected graph. Skipping...')
            
    hMap = '/pathStats'
    largest = max([len(algopaths.keys()) for algo,algopaths in collection.items()])
    allpathsizes = set()
    for algo,algopaths in collection.items():
        allpathsizes.update(set(algopaths.keys()))
    
    for algo, algopaths in collection.items():
        notpresent = allpathsizes.difference(set(algopaths.keys()))
        collection[algo].update({np:0 for np in notpresent})
                   
    dataDF = pd.DataFrame(collection)
    dataDF = dataDF.T
    dataDF.to_csv(outDir+hMap+'.csv')

   
    
    
def getNetProp(inGraph):
    '''
    Function to compute properties
    of a given network.
    '''

    # number of weakly connected components in 
    # reference network
    numCC = len(list(nx.weakly_connected_components(inGraph)))
    
    # number of feedback loop 
    # in reference network
    allCyc = nx.simple_cycles(inGraph)
    cycSet = set()
    for cyc in allCyc:
        if len(cyc) == 3:
            cycSet.add(frozenset(cyc))
    
    numFB = len(cycSet)
    
    # number of feedfwd loops
    # in reference network
    allPaths = []
    allPathsSet = set()   
    for u,v in inGraph.edges():
        allPaths = nx.all_simple_paths(inGraph, u, v, cutoff=2)
        for p in allPaths:
            if len(p) >  2:
                allPathsSet.add(frozenset(p))
                
    numFF= len(allPathsSet)
    
    
    # number of mutual interactions
    numMI = 0.0
    for u,v in inGraph.edges():
        if (v,u) in inGraph.edges():
            numMI += 0.5

    return numCC, numFB, numFF, numMI

def getEdgeHistogram(inGraph, refGraph):
    falsePositives = set(inGraph.edges()).difference(refGraph.edges())
    edgeHistogramCounts = {0:0}
    
    for fe in falsePositives:
        u,v = fe
        try:
            path = nx.dijkstra_path(refGraph,u,v)
            pathlength = len(path) -1 
            if pathlength in edgeHistogramCounts.keys():
                edgeHistogramCounts[pathlength] +=1
            else:
                edgeHistogramCounts[pathlength] = 0
            
        except nx.exception.NetworkXNoPath:
            edgeHistogramCounts[0] +=1
    return edgeHistogramCounts



def pathStats(inGraph, refGraph):
    """
    Only returns TP, FP, numPredictions for each networks
    """
    falsePositives = set(inGraph.edges()).difference(refGraph.edges())
    truePositives = set(inGraph.edges()).intersection(refGraph.edges())
    numPredictions = len(inGraph.edges())
    nopath = 0
    yespath = 0
    edgeCounts = {0:0,2:0,3:0,4:0,5:0}    
    for fe in falsePositives:
        u,v = fe
        try:
            path = nx.dijkstra_path(refGraph,u,v)
            pathlength = len(path) -1
            yespath +=1
            if pathlength in edgeCounts.keys():
                edgeCounts[pathlength] +=1
            
        except nx.exception.NetworkXNoPath:
            nopath +=1

    edgeCounts['numPred'] = numPredictions
    edgeCounts['numTP'] = len(truePositives)
    edgeCounts['numFP_withPath'] = yespath
    edgeCounts['numFP_noPath'] = nopath
    return edgeCounts
    