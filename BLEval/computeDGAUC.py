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

import logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def PRROC(dataDict, inputSettings, directed = True, selfEdges = False, plotFlag = False):
    logger.info("we have entered the dgauc file")
    '''
    Computes areas under the precision-recall and ROC curves
    for a given dataset for each algorithm.
    

    :param directed:   A flag to indicate whether to treat predictions as directed edges (directed = True) or undirected edges (directed = False).
    :type directed: bool
    :param selfEdges:   A flag to indicate whether to includeself-edges (selfEdges = True) or exclude self-edges (selfEdges = False) from evaluation.
    :type selfEdges: bool
    :param plotFlag:   A flag to indicate whether or not to save PR and ROC plots.
    :type plotFlag: bool
        
    :returns:
            - AUPRC: A dictionary containing AUPRC values for each algorithm
            - AUROC: A dictionary containing AUROC values for each algorithm
    '''
    
    # Read file for trueEdges
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', 
                                header = 0, index_col = None)
            
    # Initialize data dictionaries
    precisionDict = {}
    recallDict = {}
    FPRDict = {}
    TPRDict = {}
    AUPRC = {}
    AUROC = {}
    
    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    
    if directed:
        logger.info("We have entered the if loop")
        for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):

            # check if the output rankedEdges file exists
            if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
                 # Initialize Precsion

                predDF = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                            sep = '\t', header =  0, index_col = None)

                precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]] = computeScores(trueEdgesDF, predDF, directed = True, selfEdges = selfEdges)

            else:
                print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                      ' does not exist. Skipping...')
            PRName = '/PRplot'
            ROCName = '/ROCplot'
    else:
        logger.info("we have entered the else loop")
        for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):

            # check if the output rankedEdges file exists
            if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
                 # Initialize Precsion

                predDF = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                            sep = '\t', header =  0, index_col = None)

                precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]] = computeScores(trueEdgesDF, predDF, directed = False, selfEdges = selfEdges)

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
    return AUPRC, AUROC

# PREVIOUS IMPLEMENTATION 
# def computeScores(trueEdgesDF, predEdgeDF, 
#                   directed = True, selfEdges = True):
#     logger.info("we have entered the compute scores function")
#     '''        
#     Computes precision-recall and ROC curves
#     using scikit-learn for a given set of predictions in the 
#     form of a DataFrame.
    
#     :param trueEdgesDF:   A pandas dataframe containing the true classes.The indices of this dataframe are all possible edgesin a graph formed using the genes in the given dataset. This dataframe only has one column to indicate the classlabel of an edge. If an edge is present in the reference network, it gets a class label of 1, else 0.
#     :type trueEdgesDF: DataFrame
        
#     :param predEdgeDF:   A pandas dataframe containing the edge ranks from the prediced network. The indices of this dataframe are all possible edges.This dataframe only has one column to indicate the edge weightsin the predicted network. Higher the weight, higher the edge confidence.
#     :type predEdgeDF: DataFrame
    
#     :param directed:   A flag to indicate whether to treat predictionsas directed edges (directed = True) or undirected edges (directed = False).
#     :type directed: bool
        
#     :param selfEdges:   A flag to indicate whether to includeself-edges (selfEdges = True) or exclude self-edges (selfEdges = False) from evaluation.
#     :type selfEdges: bool
        
#     :returns:
#             - prec: A list of precision values (for PR plot)
#             - recall: A list of precision values (for PR plot)
#             - fpr: A list of false positive rates (for ROC plot)
#             - tpr: A list of true positive rates (for ROC plot)
#             - AUPRC: Area under the precision-recall curve
#             - AUROC: Area under the ROC curve
#     '''

#     if directed:        
#         logger.debug("we have entered the if statement in computescores()")
#         # Initialize dictionaries with all 
#         # possible edges
#         if selfEdges:
#             possibleEdges = list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
#                                          repeat = 2))
#         else:
#             possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
#                                          r = 2))
        
#         TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
#         PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        
#         # Compute TrueEdgeDict Dictionary
#         # 1 if edge is present in the ground-truth
#         # 0 if edge is not present in the ground-truth
        
#         loop_counter = 0
#         logger.info("loop inside the computescores() starting")
#         for key in TrueEdgeDict.keys():
#             logger.info(f"this is the {key} that we are working with")
#             if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
#                    (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
#                     logger.info(f"splitting on this {key} performed")
#                     TrueEdgeDict[key] = 1
#             loop_counter += 1
#             logger.info(f"loop count == {loop_counter}")        
                    
#         logger.info("loop exited")
                
#         logger.info("second loop started")        
#         for key in PredEdgeDict.keys():
#             subDF = predEdgeDF.loc[(predEdgeDF['Gene1'] == key.split('|')[0]) &
#                                (predEdgeDF['Gene2'] == key.split('|')[1])]
#             if len(subDF)>0:
#                 PredEdgeDict[key] = np.abs(subDF.EdgeWeight.values[0])
#         logger.info("second loop ended")
#     # if undirected
#     else:
        
#         # Initialize dictionaries with all 
#         # possible edges
#         if selfEdges:
#             possibleEdges = list(combinations_with_replacement(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
#                                                                r = 2))
#         else:
#             possibleEdges = list(combinations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
#                                                                r = 2))
#         TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
#         PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}

#         # Compute TrueEdgeDict Dictionary
#         # 1 if edge is present in the ground-truth
#         # 0 if edge is not present in the ground-truth

#         for key in TrueEdgeDict.keys():
#             if len(trueEdgesDF.loc[((trueEdgesDF['Gene1'] == key.split('|')[0]) &
#                            (trueEdgesDF['Gene2'] == key.split('|')[1])) |
#                               ((trueEdgesDF['Gene2'] == key.split('|')[0]) &
#                            (trueEdgesDF['Gene1'] == key.split('|')[1]))]) > 0:
#                 TrueEdgeDict[key] = 1  

#         # Compute PredEdgeDict Dictionary
#         # from predEdgeDF

#         for key in PredEdgeDict.keys():
#             subDF = predEdgeDF.loc[((predEdgeDF['Gene1'] == key.split('|')[0]) &
#                                (predEdgeDF['Gene2'] == key.split('|')[1])) |
#                               ((predEdgeDF['Gene2'] == key.split('|')[0]) &
#                                (predEdgeDF['Gene1'] == key.split('|')[1]))]
#             if len(subDF)>0:
#                 PredEdgeDict[key] = max(np.abs(subDF.EdgeWeight.values))

                
                
#     # Combine into one dataframe
#     # to pass it to sklearn
#     logger.info("final dataframe formed")
#     outDF = pd.DataFrame([TrueEdgeDict,PredEdgeDict]).T
#     logger.info("columns of trueEdges and PredEdges being formed")
#     outDF.columns = ['TrueEdges','PredEdges']
#     logger.info("final PRROC import taking place")
#     prroc = importr('PRROC')
#     prCurve = prroc.pr_curve(scores_class0 = FloatVector(list(outDF['PredEdges'].values)), 
#               weights_class0 = FloatVector(list(outDF['TrueEdges'].values)))

#     fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
#                                      y_score=outDF['PredEdges'], pos_label=1)

#     prec, recall, thresholds = precision_recall_curve(y_true=outDF['TrueEdges'],
#                                                       probas_pred=outDF['PredEdges'], pos_label=1)
#     logger.info("final before return")
    
#     return prec, recall, fpr, tpr, prCurve[2][0], auc(fpr, tpr)


# NEW IMPLEMENTATION


def computeScores(trueEdgesDF, predEdgeDF, 
                  directed = True, selfEdges = True):
    logger.info("we have entered the compute scores function")
    
    # Get unique genes
    unique_genes = np.unique(trueEdgesDF.loc[:, ['Gene1', 'Gene2']])
    
    # Convert DataFrames to dictionaries for faster lookup
    true_edges = set(map(tuple, trueEdgesDF[['Gene1', 'Gene2']].values))
    pred_edges = dict(zip(map(tuple, predEdgeDF[['Gene1', 'Gene2']].values), 
                          predEdgeDF['EdgeWeight']))

    TrueEdgeDict = {}
    PredEdgeDict = {}

    logger.info("Starting edge evaluation")

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

    logger.info("Edge evaluation completed")

    # Combine into one dataframe
    outDF = pd.DataFrame({'TrueEdges': TrueEdgeDict, 'PredEdges': PredEdgeDict})

    logger.info("final PRROC import taking place")
    prroc = importr('PRROC')
    prCurve = prroc.pr_curve(scores_class0 = FloatVector(list(outDF['PredEdges'].values)), 
              weights_class0 = FloatVector(list(outDF['TrueEdges'].values)))

    fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
                                     y_score=outDF['PredEdges'], pos_label=1)

    prec, recall, thresholds = precision_recall_curve(y_true=outDF['TrueEdges'],
                                                      probas_pred=outDF['PredEdges'], pos_label=1)
    logger.info("final before return")
    
    return prec, recall, fpr, tpr, prCurve[2][0], auc(fpr, tpr)