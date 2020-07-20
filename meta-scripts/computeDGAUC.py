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
from rpy2.robjects.packages import importr
from rpy2.robjects import FloatVector


def PRROC(dataDict, inputSettings, directed = True, TFEdges = False, plotFlag = False):
    '''
    Computes areas under the precision-recall and ROC curves
    for a given dataset for each algorithm.
    
    Parameters
    -----------
        directed: bool
            A flag to indicate whether to treat predictions
            as directed edges (directed = True) or 
            undirected edges (directed = False).
            
        TFEdges: bool
            A flag to indicate whether to include
            self-edges (TFEdges = True) or 
            include only TF-gene genes (TFEdges = True) for evaluation.
        
        plotFlag: bool
            A flag to indicate whether or not to save PR and ROC plots.
            
    :returns:
            - AUPRC: A dictionary containing AUPRC values for each algorithm
            - AUROC: A dictionary containing AUROC values for each algorithm
    '''
    
    # Read file for trueEdges
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', 
                                header = 0, index_col = None)
    trueEdgesDF = trueEdgesDF.loc[(trueEdgesDF['Gene1'] != trueEdgesDF['Gene2'])]
    trueEdgesDF.drop_duplicates(keep = 'first', inplace=True)
    trueEdgesDF.reset_index(drop=True, inplace=True)

    # Initialize data dictionaries
    precisionDict = {}
    recallDict = {}
    FPRDict = {}
    TPRDict = {}
    AUPRC = {}
    AUROC = {}
    
    trueEdgeCnt = {}
    possEdgeCnt = {}
    TFCnt = {}
    nodeCnt = {}
    
    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    
    if directed:
        for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):

            # check if the output rankedEdges file exists
            if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
                 # Initialize Precsion

                predDF = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                            sep = '\t', header =  0, index_col = None)
                predDF.fillna(0, inplace = True)
                precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]], trueEdgeCnt[algo[0]], possEdgeCnt[algo[0]], TFCnt[algo[0]], nodeCnt[algo[0]]  = computeScores(trueEdgesDF, predDF, directed = True, TFEdges = TFEdges)

            else:
                print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                      ' does not exist. Skipping...')
            PRName = '/PRplot'
            ROCName = '/ROCplot'
    else:
        for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):

            # check if the output rankedEdges file exists
            if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
                 # Initialize Precsion

                predDF = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                            sep = '\t', header =  0, index_col = None)

                precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]] = computeScores(trueEdgesDF, predDF, directed = False, TFEdges = TFEdges)

            else:
                print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' does not exist. Skipping...')
            
            PRName = '/uPRplot'
            ROCName = '/uROCplot'
    if (plotFlag):
         ## Make PR curves
        legendList = []
        for key in recallDict.keys():
            sns.lineplot(recallDict[key],precisionDict[key], ci=None)
            legendList.append(key + ' (AUPRC = ' + str("%.2f" % (AUPRC[key]))+')')
        plt.xlim(0,1)    
        plt.ylim(0,1)
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.legend(legendList) 
        plt.savefig(outDir+PRName+'.pdf')
        plt.savefig(outDir+PRName+'.png')
        plt.clf()

        ## Make ROC curves
        legendList = []
        for key in recallDict.keys():
            sns.lineplot(FPRDict[key],TPRDict[key], ci=None)
            legendList.append(key + ' (AUROC = ' + str("%.2f" % (AUROC[key]))+')')

        plt.plot([0, 1], [0, 1], linewidth = 1.5, color = 'k', linestyle = '--')

        plt.xlim(0,1)    
        plt.ylim(0,1)
        plt.xlabel('FPR')
        plt.ylabel('TPR')
        plt.legend(legendList) 
        plt.savefig(outDir+ROCName+'.pdf')
        plt.savefig(outDir+ROCName+'.png')
        plt.clf()
    return AUPRC, AUROC, trueEdgeCnt, possEdgeCnt, TFCnt, nodeCnt

def computeScores(trueEdgesDF, predEdgeDF, 
                  directed = True, TFEdges = True):
    '''        
    Computes precision-recall and ROC curves
    using scikit-learn for a given set of predictions in the 
    form of a DataFrame.
    
    Parameters
    -----------
        trueEdgesDF: DataFrame
            A pandas dataframe containing the true classes.
            The indices of this dataframe are all possible edges
            in a graph formed using the genes in the given dataset. 
            This dataframe only has one column to indicate the class
            label of an edge. If an edge is present in the reference
            network, it gets a class label of 1, else 0.
            
        predEdgeDF: DataFrame
            A pandas dataframe containing the edge ranks from the prediced 
            network. The indices of this dataframe are all possible edges.
            This dataframe only has one column to indicate the edge weights
            in the predicted network. Higher the weight, higher the 
            edge confidence.
        
        directed: bool
            A flag to indicate whether to treat predictions
            as directed edges (directed = True) or 
            undirected edges (directed = False).
            
        TFEdges: bool
            A flag to indicate whether to include
            self-edges (TFEdges = True) or 
            exclude self-edges (TFEdges = False) from evaluation.
            
    :returns:
            - prec: A list of precision values (for PR plot)
            - recall: A list of precision values (for PR plot)
            - fpr: A list of false positive rates (for ROC plot)
            - tpr: A list of true positive rates (for ROC plot)
            - AUPRC: Area under the precision-recall curve
            - AUROC: Area under the ROC curve
    '''

    if directed:        
        # Initialize dictionaries with all 
        # possible edges
        if TFEdges:
            
            # Get a list of all possible TF to gene interactions 
            uniqueNodes = np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']])
            possibleEdges_TF = set(product(set(trueEdgesDF.Gene1),set(uniqueNodes)))

            # Get a list of all possible interactions 
            possibleEdges_noSelf = set(permutations(uniqueNodes, r = 2))
            
            # Find intersection of above lists to ignore self edges
            # TODO: is there a better way of doing this?
            possibleEdges = possibleEdges_TF.intersection(possibleEdges_noSelf)
            #print("\n Total:",len(set(trueEdgesDF.Gene1)), "possible TFs" )
            #print("\n Total:",len(set(trueEdgesDF.Gene2)), "possible Target Genes" )

            #print("\n Total:",len(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']])), "possible nodes" )
            #print("\n Total:",len(possibleEdges), "possible edges")
            
        else:
            uniqueNodes = np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']])

            possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                         r = 2))

        TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}

        trueEdges = trueEdgesDF['Gene1'] + "|" + trueEdgesDF['Gene2']
        trueEdges = trueEdges[trueEdges.isin(TrueEdgeDict)]
        predEdgeDF['Edges'] = predEdgeDF['Gene1'] + "|" + predEdgeDF['Gene2']
        # limit the predicted edges to the genes that are in the ground truth
        predEdgeDF = predEdgeDF[predEdgeDF['Edges'].isin(TrueEdgeDict)]

        # Compute TrueEdgeDict Dictionary
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth
        trueEdgeCnt = 0
        for edge in trueEdges:
            TrueEdgeDict[edge] = 1
            trueEdgeCnt += 1
        for edge, score in zip(predEdgeDF['Edges'], predEdgeDF['EdgeWeight']):
            PredEdgeDict[edge] = np.abs(score)

    # if undirected
    else:
        
        # Initialize dictionaries with all 
        # possible edges
        if TFEdges:
            possibleEdges = list(product(set(trueEdgesDF.Gene1),set(trueEdgesDF.Gene2)))

        else:
            possibleEdges = list(combinations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                                               r = 2))
        TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}

        trueEdges = trueEdgesDF['Gene1'] + "|" + trueEdgesDF['Gene2']
        trueEdges = trueEdges[trueEdges.isin(TrueEdgeDict)]

        # Compute PredEdgeDict Dictionary
        # from predEdgeDF

        for key in PredEdgeDict.keys():
            subDF = predEdgeDF.loc[((predEdgeDF['Gene1'] == key.split('|')[0]) &
                               (predEdgeDF['Gene2'] == key.split('|')[1])) |
                              ((predEdgeDF['Gene2'] == key.split('|')[0]) &
                               (predEdgeDF['Gene1'] == key.split('|')[1]))]
            if len(subDF)>0:
                PredEdgeDict[key] = max(np.abs(subDF.EdgeWeight.values))

                
                
    # Combine into one dataframe
    # to pass it to sklearn
    outDF = pd.DataFrame([TrueEdgeDict,PredEdgeDict]).T
    outDF.columns = ['TrueEdges','PredEdges']
    prroc = importr('PRROC')

    prCurve = prroc.pr_curve(scores_class0 = FloatVector(list(outDF['PredEdges'].values)), 
              weights_class0 = FloatVector(list(outDF['TrueEdges'].values)))

    fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
                                     y_score=outDF['PredEdges'], pos_label=1)

    prec, recall, thresholds = precision_recall_curve(y_true=outDF['TrueEdges'],
                                                      probas_pred=outDF['PredEdges'], pos_label=1)

    #print("\nUsing PRROC: ", prCurve[2][0]/(trueEdgeCnt/len(possibleEdges)))
    print("\nEdges considered ",len(trueEdges))
    return prec, recall, fpr, tpr, prCurve[2][0], auc(fpr, tpr), len(trueEdges), len(possibleEdges), len(set(trueEdgesDF.Gene1)), len(uniqueNodes)
