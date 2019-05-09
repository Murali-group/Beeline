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

def motifAnalysis(dataDict, inputSettings):
    '''
    Computes "directed","feed-forward", 
    "cascade", and "mutual" motifs.
    '''
    
    # Read file for trueEdges
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', 
                                header = 0, index_col = None)
            
    
    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    dataDict = {}
    dataDict['Directionality'] = {}
    dataDict['FFL'] = {}
    dataDict['Cascade'] = {}
    dataDict['FanIn'] = {}
    dataDict['FanOut'] = {}

#    for algo in tqdm(inputSettings.algorithms, 
#                     total = len(inputSettings.algorithms), unit = " Algorithms"):
    for algo in inputSettings.algorithms:
        # check if the output rankedEdges file exists
        if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
             # Initialize Precsion

            predDF = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                        sep = '\t', header =  0, index_col = None)
            # min-max normalize edge weights to have a range of 0-1
            predDF['EdgeWeight'] = (predDF['EdgeWeight']-predDF['EdgeWeight'].min())/(predDF['EdgeWeight'].max()-predDF['EdgeWeight'].min())
            dataDict['Directionality'][algo[0]], dataDict['FFL'][algo[0]], dataDict['Cascade'][algo[0]], dataDict['FanIn'][algo[0]], dataDict['FanOut'][algo[0]]  = getMotifScore(trueEdgesDF, predDF)

        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' does not exist. Skipping...')
        hMap = '/motifRes'
        
     ## Make heatmap
    palette = sns.diverging_palette(10, 220, sep=1, as_cmap = True)

    dataDF = pd.DataFrame(dataDict)
    #dataDF = (dataDF-dataDF.min())/(dataDF.max()-dataDF.min())

    dataDF.to_csv(outDir+hMap+'.csv')
    legendList = []
    maxVal = max(dataDF.max().max(),abs(dataDF.min().min()))
    sns.clustermap(dataDF, cmap = palette, vmax= 0.3, vmin=-0.3)
    plt.savefig(outDir+hMap+'.pdf')
    plt.savefig(outDir+hMap+'.png')
    plt.clf()

def getMotifScore(trueEdgesDF, predEdgeDF):
    '''
    Function to compute a score for all motifs.
    '''

    possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                 r = 2))        
    EdgeDict = {'|'.join(p):0 for p in possibleEdges}

    refGraph = nx.DiGraph()
    # Compute TrueEdgeDict Dictionary
    # 1 if edge is present in the ground-truth
    # 0 if edge is not present in the ground-truth

    for key in EdgeDict.keys():
        u = key.split('|')[0]
        v = key.split('|')[1]
        if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == u) &
               (trueEdgesDF['Gene2'] == v)])>0:
                refGraph.add_edge(u,v)

    predGraph = nx.DiGraph()

    for key in EdgeDict.keys():
        u = key.split('|')[0]
        v = key.split('|')[1]
        subDF = predEdgeDF.loc[(predEdgeDF['Gene1'] == key.split('|')[0]) &
                               (predEdgeDF['Gene2'] == key.split('|')[1])]
        if len(subDF)>0:
            predGraph.add_edge(u,v, weight = np.abs(subDF.EdgeWeight.values[0]))
        else:
            predGraph.add_edge(u,v, weight = 0)

    # Initialize motif scores
    r_direction = 0
    r_feedfor = 0
    r_casc = 0
    r_mutual = 0

    # Directionlity motif
    # Compute bias for finding correct
    # directionality of an edge
    rm = 0.0
    rm_bar = 0.0
    for u,v in refGraph.edges():
        data = predGraph.get_edge_data(u,v)
        rm += data['weight']
        data = predGraph.get_edge_data(v,u)
        rm_bar += data['weight']
    r_direction = (rm/len(refGraph.edges()) - rm_bar/len(refGraph.edges()))

    # Feed-forward loop motif
    # Compute bias for finding correct
    # directionality of an edge
    rm = 0.0
    rm_bar = 0.0
    for u,v in refGraph.edges():
        allPaths = []
        allPathsList = []
        allPaths = nx.all_simple_paths(refGraph, u, v, cutoff=2)
        allPathsList = [p for p in allPaths]
        # first path in allPathsList is u,v
        # if paths of at most length 2 exist other than u,v:
        if len(allPathsList) > 1:
            data = predGraph.get_edge_data(u,v)
            rm += data['weight']
        else:
            data = predGraph.get_edge_data(u,v)
            rm_bar += data['weight']
    r_feedfor = (rm/len(refGraph.edges()) - rm_bar/len(refGraph.edges()))

    
    # Bias to predict false positives A->B
    #where an indirect path A->C->B
    #exists (vs. false positives where no
    #indirect path exists)
    rm = 0.0
    rm_bar = 0.0
    refGraphComp = nx.complement(refGraph)
    for u,v in refGraphComp.edges():
        allPaths = []
        allPathsList = []
        allPaths = nx.all_simple_paths(refGraph, u, v, cutoff=2)
        allPathsList = [p for p in allPaths]
        # first path in allPathsList is u,v
        # if paths of at most length 2 exist other than u,v:
        if len(allPathsList) > 0:
            data = predGraph.get_edge_data(u,v)
            rm += data['weight']
        else:
            data = predGraph.get_edge_data(u,v)
            rm_bar += data['weight']
    r_casc = (rm/len(refGraph.edges()) - rm_bar/len(refGraph.edges()))
    
    
    # Fan-in motif
    # Compute bias for finding correct
    # directionality of an edge
    rm = 0.0
    rm_bar = 0.0
    for u,v in refGraph.edges():
        prede = [p for p in refGraph.predecessors(v)]
        if len(prede)>1:
            data = predGraph.get_edge_data(u,v)
            rm += data['weight']
        else:
            data = predGraph.get_edge_data(u,v)
            rm_bar += data['weight']
    r_fanin = (rm/len(refGraph.edges()) - rm_bar/len(refGraph.edges()))
    
    
    # Fan-out motif
    # Compute bias for finding correct
    # directionality of an edge
    rm = 0.0
    rm_bar = 0.0
    for u,v in refGraph.edges():
        succe = [p for p in refGraph.successors(u)]
        if len(succe)>1:
            data = predGraph.get_edge_data(u,v)
            rm += data['weight']
        else:
            data = predGraph.get_edge_data(u,v)
            rm_bar += data['weight']
    r_fanout = (rm/len(refGraph.edges()) - rm_bar/len(refGraph.edges()))
    
    return r_direction, r_feedfor, r_casc, r_fanin, r_fanout
