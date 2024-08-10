import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, roc_auc_score
from sklearn.metrics import f1_score #, accuracy_score
from itertools import product, permutations, combinations, combinations_with_replacement
from tqdm import tqdm
from rpy2.robjects.packages import importr
from rpy2.robjects import FloatVector


def signedPRROC(dataDict, inputSettings, directed = True, selfEdges = False, plotFlag = False):
    '''
    Computes areas under the precision-recall and ROC curves of activation edges and inhibitory edges 
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
            - AP: A dictionary containing AP values for each algorithm
            - F1: A dictionary containing F1 values for each algorithm
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
    AP = {}
    F1 = {}
    
    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    
    # Obtation predicted dataframe for each algorithm
    for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):
        # check if the output rankedEdges file exists
        if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
            
            # Initialize Precsion
            predDF = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                        sep = '\t', header =  0, index_col = None)

            precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]], AP[algo[0]], F1[algo[0]] = signedComputeScores(trueEdgesDF, predDF, directed = directed, selfEdges = selfEdges)
        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                    ' does not exist. Skipping...')
       
        PRName = '/PRplot'
        ROCName = '/ROCplot'
        APName = '/APplot'
                

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

        ## Make AP curves
        legendList = []
        for key in recallDict.keys():
            sns.lineplot(recallDict[key],precisionDict[key], ci=None)
            legendList.append(key + ' (AP = ' + str("%.2f" % (AP[key]))+')')
        plt.xlim(0,1)    
        plt.ylim(0,1)
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.legend(legendList) 
        plt.savefig(outDir+APName+'.pdf')
        plt.savefig(outDir+APName+'.png')
        plt.clf()

    return AUPRC, AUROC, AP, F1


def signedComputeScores(trueEdgesDF, predEdgeDF, 
                  directed = True, selfEdges = True):
    '''        
    Computes precision-recall and ROC curves
    using scikit-learn for a given set of predictions in the 
    form of a DataFrame.
    
    :param trueEdgesDF:   A pandas dataframe containing the true classes.The indices of this dataframe are all possible edges in a graph formed using the genes in the given dataset. This dataframe only has one column to indicate the class label of an edge. If an edge is present in the reference network, it gets a class label of 1, else 0.
    :type trueEdgesDF: DataFrame
        
    :param predEdgeDF:   A pandas dataframe containing the edge ranks from the predicted network. The indices of this dataframe are all possible edges. This dataframe only has one column to indicate the edge weights in the predicted network. Higher the weight, higher the edge confidence.
    :type predEdgeDF: DataFrame
    
    :param directed:   A flag to indicate whether to treat predictionsas directed edges (directed = True) or undirected edges (directed = False).
    :type directed: bool
        
    :param selfEdges:   A flag to indicate whether to include self-edges (selfEdges = True) or exclude self-edges (selfEdges = False) from evaluation.
    :type selfEdges: bool
        
    :returns:
            - prec: A list of precision values (for PR plot)
            - recall: A list of precision values (for PR plot)
            - fpr: A list of false positive rates (for ROC plot)
            - tpr: A list of true positive rates (for ROC plot)
            - AUPRC: Area under the precision-recall curve
            - AUROC: Area under the ROC curve
            - AP: Average precision
            - F1: F1 score
    '''

    # Create lists for possible directed and undirected edges
    if directed:
        if selfEdges:
            possibleEdges = list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                            repeat = 2))
        else:
            # permutations of gene pairs
            possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                            r = 2))
    else:
        if selfEdges:
            possibleEdges = list(combinations_with_replacement(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                                                r = 2))
        else:
            # combination of gene pairs
            possibleEdges = list(combinations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]), r = 2))
    
    # Determine if the prediction is signed using nonnegative edge weights
    is_pred_signed = np.count_nonzero(predEdgeDF.EdgeWeight >= 0) != len(predEdgeDF)                                                             
    outDFAll = {'+':{},'-':{}}

    # Consider signs
    for sgn in ['+', "-"]:
        # Initialize dictionaries with all possible edges
        # Obtain the dictionary of gene pairs for true edges and predicted edges and assign 0 for all pairs
        TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}
    
        # Store edges with different sign
        ignoredEdges = set()

        # Compute TrueEdgeDict Dictionary
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth
        for edge in trueEdgesDF.itertuples():
            if edge.Gene1 == edge.Gene2:
                continue
            
            if edge.Type == sgn:
                if "|".join((edge.Gene1, edge.Gene2)) in TrueEdgeDict:
                    TrueEdgeDict["|".join((edge.Gene1, edge.Gene2))] = 1
                
                if not directed:
                    if "|".join((edge.Gene2, edge.Gene1)) in TrueEdgeDict:
                        TrueEdgeDict["|".join((edge.Gene2, edge.Gene1))] = 1

            else:
                # Ignored edges not in ground-truth or with diffrent sign
                ignoredEdges.add("|".join((edge.Gene1, edge.Gene2)))

        # Compute PredEdgeDict Dictionary
        for edge in predEdgeDF.itertuples():
            if edge.Gene1 == edge.Gene2:
                continue

            if is_pred_signed:
                # Determine signs based on predicted edge weights
                edge_sign = "+" if edge.EdgeWeight >= 0 else "-"
                # Absolute value of predicted edge weight if edge is with the same sign of interest
                if edge_sign == sgn:
                    if "|".join((edge.Gene1, edge.Gene2)) in PredEdgeDict:
                        PredEdgeDict["|".join((edge.Gene1, edge.Gene2))] = np.abs(edge.EdgeWeight)

            else:
                if "|".join((edge.Gene1, edge.Gene2)) in ignoredEdges:
                    continue
                # Assign absoulute values of edge weights to edges if the predicted edges are not signed
                PredEdgeDict["|".join((edge.Gene1, edge.Gene2))] = np.abs(edge.EdgeWeight)
           
        outDF = pd.DataFrame([TrueEdgeDict,PredEdgeDict]).T
        outDF.columns = ['TrueEdges','PredEdges']
        outDFAll[sgn] = outDF

    # Combine into one dataframe and pass to sklearn
    prroc = importr('PRROC')
    precAll = {'+':{},'-':{}}
    recallAll = {'+':{},'-':{}}
    fprAll = {'+':{},'-':{}}
    tprAll = {'+':{},'-':{}}
    auprcAll = {'+':{},'-':{}}
    aurocAll = {'+':{},'-':{}}
    apAll = {'+':{}, '-':{}}
    f1All = {'+':{}, '-':{}}
    
    # Compute AUC by sign
    for sgn in ['+','-']:
        # precision, recall
        prec, recall, thresholds = precision_recall_curve(y_true=outDFAll[sgn]['TrueEdges'],
                                                          probas_pred=outDFAll[sgn]['PredEdges'], 
                                                          pos_label=1)
        # FPR, TPR
        fpr, tpr, thresholds = roc_curve(y_true=outDFAll[sgn]['TrueEdges'],
                                         y_score=outDFAll[sgn]['PredEdges'], 
                                         pos_label=1)
        # AUPRC
        auprc = prroc.pr_curve(scores_class0 = FloatVector(list(outDFAll[sgn]['PredEdges'].values)), 
                               weights_class0 = FloatVector(list(outDFAll[sgn]['TrueEdges'].values)))
        # AUROC
        auroc = roc_auc_score(y_true=outDFAll[sgn]['TrueEdges'],
                              y_score=outDFAll[sgn]['PredEdges'])
        
        # AP
        ap = average_precision_score(y_true=outDFAll[sgn]['TrueEdges'], 
                                     y_score=outDFAll[sgn]['PredEdges'], 
                                     pos_label=1)
        
        # We use the weighted average on the data.
        # Please refer to the selection of avaerage options: https://towardsdatascience.com/micro-macro-weighted-averages-of-f1-score-clearly-explained-b603420b292f
        f1 = f1_score(y_true=outDF['TrueEdges'], 
                    y_pred=outDF['PredEdges'].round(), 
                    pos_label=1, average='weighted')
        
        # # The accuracy is the same as f1 with micro average
        # accuracy = accuracy_score(y_true=outDF['TrueEdges'], 
        #                           y_pred=outDF['PredEdges'].round())
                                
        precAll[sgn] = prec
        recallAll[sgn] = recall
        fprAll[sgn] = fpr
        tprAll[sgn] = tpr
        auprcAll[sgn] = auprc[2][0]
        aurocAll[sgn] = auroc
        apAll[sgn] = ap
        f1All[sgn] = f1

    return precAll, recallAll, fprAll, tprAll, auprcAll, aurocAll, apAll, f1All
