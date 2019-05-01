import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from itertools import product, combinations_with_replacement

def ComputeEnsemble(dataDict, inputSettings, directed = True):
    '''
    Computes ensemble network
    for the given set of algorithms
    '''
    
    # Read file for trueEdges
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', header = 0, index_col = None)
    
    TrueEdgeDict = {'|'.join(p):0 for p in list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),repeat =2))}
    for key in TrueEdgeDict.keys():
        if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
            TrueEdgeDict[key] = 1
            def ComputeEnsemble(dataDict, inputSettings):
    '''
    Computes ensemble network
    for the given set of algorithms
    '''
    
    # Read file for trueEdges
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', header = 0, index_col = None)
    
    TrueEdgeDict = {'|'.join(p):0 for p in list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),repeat =2))}
    for key in TrueEdgeDict.keys():
        if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
            TrueEdgeDict[key] = 1
            
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    
    
    Pred = {}
    for algo in inputSettings.algorithms:
        if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
            Pred[algo[0]] = {'|'.join(p):0 for p in list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                                                          repeat =2))}
            predDF =  pd.read_csv(outDir+'/'+algo[0]+'/rankedEdges.csv',sep='\t')
            for key in Pred[algo[0]].keys():
                subDF = predDF.loc[(predDF['Gene1'] == key.split('|')[0]) &
                                   (predDF['Gene2'] == key.split('|')[1])]
                if len(subDF)>0:
                    Pred[algo[0]][key] = np.abs(subDF.EdgeWeight.values[0])
                        
    ensembleDF = pd.DataFrame(Pred)
    print(ensembleDF.head())

    ensembleFinal = pd.DataFrame(ensembleDF.rank(axis = 'index', ascending = True).mean(axis= 'columns'))
    print(ensembleFinal.head())

    ensembleFinal.loc[:,'Gene1'] = [ix.split('|')[0] for ix in ensembleFinal.index]
    ensembleFinal.loc[:,'Gene2'] = [ix.split('|')[1] for ix in ensembleFinal.index]
    ensembleFinal.columns = ['EdgeWeight','Gene1','Gene2']
    ensembleFinal.sort_values(by='EdgeWeight',ascending = False, inplace = True)
    print("Writing Ensemble Model output to ",outDir)
    ensembleFinal[['Gene1','Gene2','EdgeWeight']].to_csv(outDir+'/rankedEdges.csv', sep = '\t', index = False)            


    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    
    
    Pred = {}
    for algo in inputSettings.algorithms:
        if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
            Pred[algo[0]] = {'|'.join(p):0 for p in list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                                                          repeat =2))}
            predDF =  pd.read_csv(outDir+'/'+algo[0]+'/rankedEdges.csv',sep='\t')
            for key in Pred[algo[0]].keys():
                subDF = predDF.loc[(predDF['Gene1'] == key.split('|')[0]) &
                                   (predDF['Gene2'] == key.split('|')[1])]
                if len(subDF)>0:
                    Pred[algo[0]][key] = np.abs(subDF.EdgeWeight.values[0])
                        
    ensembleDF = pd.DataFrame(Pred)
    print(ensembleDF.head())

    ensembleFinal = pd.DataFrame(ensembleDF.rank(axis = 'index', ascending = True).mean(axis= 'columns'))
    print(ensembleFinal.head())

    ensembleFinal.loc[:,'Gene1'] = [ix.split('|')[0] for ix in ensembleFinal.index]
    ensembleFinal.loc[:,'Gene2'] = [ix.split('|')[1] for ix in ensembleFinal.index]
    ensembleFinal.columns = ['EdgeWeight','Gene1','Gene2']
    ensembleFinal.sort_values(by='EdgeWeight',ascending = False, inplace = True)
    print("Writing Ensemble Model output to ",outDir)
    ensembleFinal[['Gene1','Gene2','EdgeWeight']].to_csv(outDir+'/rankedEdges.csv', sep = '\t', index = False)            

