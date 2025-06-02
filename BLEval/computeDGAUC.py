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


def PRROC(dataDict, inputSettings, directed=True, selfEdges=False, plotFlag=False):
    '''
    Computes areas under the precision-recall and ROC curves
    for a given dataset for each algorithm.
    
    :param directed:   A flag to indicate whether to treat predictions as directed edges (directed = True) or undirected edges (directed = False).
    :type directed: bool
    :param selfEdges:   A flag to indicate whether to include self-edges (selfEdges = True) or exclude self-edges (selfEdges = False) from evaluation.
    :type selfEdges: bool
    :param plotFlag:   A flag to indicate whether or not to save PR and ROC plots.
    :type plotFlag: bool
        
    :returns:
            - AUPRC: A dictionary containing AUPRC values for each algorithm
            - AUROC: A dictionary containing AUROC values for each algorithm
    '''
    
    if inputSettings.use_embeddings == True:     
        true_edges_path = inputSettings.get_true_edges_path(dataDict['name'])
        trueEdgesDF = pd.read_csv(true_edges_path, sep=',', header=0, index_col=None)
        outDir = Path("outputs") / inputSettings.datadir.relative_to("inputs") / dataDict['name'] / 'processed_ExpressionData'
    else:
        trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', 
                                header = 0, index_col = None)
        outDir = Path("outputs") / inputSettings.datadir.relative_to("inputs") / dataDict['name']

    
    # Initialize data dictionaries
    precisionDict = {}
    recallDict = {}
    FPRDict = {}
    TPRDict = {}
    AUPRC = {}
    AUROC = {}
    
    
    if directed:
        for algo in tqdm(inputSettings.algorithms, 
                         total=len(inputSettings.algorithms), unit=" Algorithms"):

            # check if the output rankedEdges file exists
            ranked_edges_path = outDir / algo[0] / 'rankedEdges.csv'
            if ranked_edges_path.exists():
                # Initialize Precision
                predDF = pd.read_csv(ranked_edges_path, sep='\t', header=0, index_col=None)

                precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]] = computeScores(trueEdgesDF, predDF, directed=True, selfEdges=selfEdges)

            else:
                print(f'{ranked_edges_path} does not exist. Skipping...')
            PRName = 'PRplot'
            ROCName = 'ROCplot'
    else:
        for algo in tqdm(inputSettings.algorithms, 
                         total=len(inputSettings.algorithms), unit=" Algorithms"):

            # check if the output rankedEdges file exists
            ranked_edges_path = outDir / algo[0] / 'rankedEdges.csv'
            if ranked_edges_path.exists():
                # Initialize Precision
                predDF = pd.read_csv(ranked_edges_path, sep='\t', header=0, index_col=None)

                precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]] = computeScores(trueEdgesDF, predDF, directed=False, selfEdges=selfEdges)

            else:
                print(f'{ranked_edges_path} does not exist. Skipping...')
            
            PRName = 'uPRplot'
            ROCName = 'uROCplot'

    if plotFlag:
        # Make PR curves
        legendList = []
        for key in recallDict.keys():
            sns.lineplot(recallDict[key], precisionDict[key], ci=None)
            legendList.append(f'{key} (AUPRC = {AUPRC[key]:.2f})')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.legend(legendList)
        plt.savefig(outDir / f'{PRName}.pdf')
        plt.savefig(outDir / f'{PRName}.png')
        plt.clf()

        # Make ROC curves
        legendList = []
        for key in recallDict.keys():
            sns.lineplot(FPRDict[key], TPRDict[key], ci=None)
            legendList.append(f'{key} (AUROC = {AUROC[key]:.2f})')

        plt.plot([0, 1], [0, 1], linewidth=1.5, color='k', linestyle='--')

        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xlabel('FPR')
        plt.ylabel('TPR')
        plt.legend(legendList)
        plt.savefig(outDir / f'{ROCName}.pdf')
        plt.savefig(outDir / f'{ROCName}.png')
        plt.clf()

    return AUPRC, AUROC

def computeScores(trueEdgesDF, predEdgeDF, 
                  directed = True, selfEdges = True):
    
    # Get unique genes
    unique_genes = np.unique(trueEdgesDF.loc[:, ['Gene1', 'Gene2']])
    
    # Convert DataFrames to dictionaries for faster lookup
    true_edges = set(map(tuple, trueEdgesDF[['Gene1', 'Gene2']].values))
    pred_edges = dict(zip(map(tuple, predEdgeDF[['Gene1', 'Gene2']].values), 
                          predEdgeDF['EdgeWeight']))

    TrueEdgeDict = {}
    PredEdgeDict = {}

    if directed:
        edge_generator = product(unique_genes, repeat=2) if selfEdges else permutations(unique_genes, r=2)
    else:
        edge_generator = combinations_with_replacement(unique_genes, r=2) if selfEdges else combinations(unique_genes, r=2)

    for edge in edge_generator:
        key = '|'.join(edge)
        
        # Compute TrueEdgeDict
        TrueEdgeDict[key] = int(edge in true_edges or (not directed and edge[::-1] in true_edges))
        
        # Compute PredEdgeDict
        if directed:
            PredEdgeDict[key] = abs(pred_edges.get(edge, 0))
        else:
            PredEdgeDict[key] = max(abs(pred_edges.get(edge, 0)), abs(pred_edges.get(edge[::-1], 0)))

    # Combine into one dataframe
    outDF = pd.DataFrame({'TrueEdges': TrueEdgeDict, 'PredEdges': PredEdgeDict})

    prroc = importr('PRROC')
    prCurve = prroc.pr_curve(scores_class0 = FloatVector(list(outDF['PredEdges'].values)), 
              weights_class0 = FloatVector(list(outDF['TrueEdges'].values)))

    fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
                                     y_score=outDF['PredEdges'], pos_label=1)

    prec, recall, thresholds = precision_recall_curve(y_true=outDF['TrueEdges'],
                                                      probas_pred=outDF['PredEdges'], pos_label=1)
    
    return prec, recall, fpr, tpr, prCurve[2][0], auc(fpr, tpr)
