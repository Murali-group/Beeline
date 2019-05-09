import pandas as pd
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

def outputAnalysis(dataDict, inputSettings):
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

    refCC, refFB, refFF, refMI = getNetProp(refGraph)

    
    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    dataDict = {}
    dataDict['Conn. Comp'] = {}
    dataDict['FFL'] = {}
    dataDict['FBL'] = {}
    dataDict['Mutual'] = {}

    for algo in tqdm(inputSettings.algorithms, 
                     total = len(inputSettings.algorithms), unit = " Algorithms"):
        
        # check if the output rankedEdges file exists
        if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
             # Initialize Precsion

            predDF = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                        sep = '\t', header =  0, index_col = None)
            
            predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
            edgeWeightTopk = predDF.iloc[len(refGraph.edges())-1].EdgeWeight
            newDF = predDF.loc[(predDF['EdgeWeight'] >= edgeWeightTopk)]
                
            predGraph = nx.DiGraph()
            

            for key in EdgeDict.keys():
                u = key.split('|')[0]
                v = key.split('|')[1]
                if len(newDF.loc[(newDF['Gene1'] == u) &
                       (newDF['Gene2'] == v)])>0:
                        predGraph.add_edge(u,v)
            
            dataDict['Conn. Comp'][algo[0]], dataDict['FBL'][algo[0]], dataDict['FFL'][algo[0]], dataDict['Mutual'][algo[0]] = getNetProp(predGraph)
            if algo[0] == 'PPCOR' or algo[0] == 'PIDC':
                dataDict['FBL'][algo[0]] = refFB
                dataDict['FFL'][algo[0]] = refFF
                dataDict['Mutual'][algo[0]] = refMI
        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' does not exist. Skipping...')
        hMap = '/topK'
        
     ## Make heatmap
    palette = sns.diverging_palette(10, 220, sep=10, as_cmap = True)

    dataDF = pd.DataFrame(dataDict)
    dataDF['Conn. Comp'] = dataDF['Conn. Comp'] - refCC
    dataDF['FFL'] = dataDF['FFL'] - refFF
    dataDF['FBL'] = dataDF['FBL'] - refFB
    dataDF['Mutual'] = dataDF['Mutual'] - refMI

    #dataDF = (dataDF-dataDF.min())/(dataDF.max()-dataDF.min())

    dataDF.to_csv(outDir+hMap+'.csv')
    #maxVal = max(dataDF.max().max(),abs(dataDF.min().min()))
    sns.heatmap(dataDF, cmap = palette, center=0)
    plt.savefig(outDir+hMap+'.pdf')
    plt.savefig(outDir+hMap+'.png')
    plt.clf()
    sys.exit()
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
            if len(p) > 2:
                allPathsSet.add(frozenset(p))
                
    numFF= len(allPathsSet)
    
    
    # number of mutual interactions
    numMI = 0.0
    for u,v in inGraph.edges():
        if (v,u) in inGraph.edges():
            numMI += 0.5

    return numCC, numFB, numFF, numMI