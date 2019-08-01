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

def PRROC(dataDict, inputSettings, directed = True, selfEdges = False,
            plotFlag = False, userReferenceNetworkFile=None,
            tfsFile=None, onlyEdgesFromTFs=False):
    '''
    Computes areas under the precision-recall and ROC curves
    for a given dataset for each algorithm.

    Parameters
    -----------
        directed: bool
            A flag to indicate whether to treat predictions
            as directed edges (directed = True) or
            undirected edges (directed = False).

        selfEdges: bool
            A flag to indicate whether to include
            self-edges (selfEdges = True) or
            exclude self-edges (selfEdges = False) from evaluation.

        plotFlag: bool
            A flag to indicate whether or not to save PR and ROC plots.

        userReferenceNetworkFile: str
            The path to a file that specifiy reference network to be used for
            AUC calculations. Default is None. If the value is not None, the
            function will overide the reference network with the given file.
            The file should be comma separated and have following node column
            names: `Gene1` and `Gene2`.

        tfsFile: str
            The path to a file that specifiy a list of transcription factors in
            the reference network. Default is None.

        onlyEdgesFromTFs: bool
            A flag to indicate whether to ignore edges from non-transcription
            factors or not. If onlyEdgesFromTFs=True, the function will try to
            fetch list of transcription factors from the file specified by
            `tfsFile` attribute and then edges from nodes other than transcription
            factors in the reference network will be ignored and considered a
            negative. Default is False.

    :returns:
            - AUPRC: A dictionary containing AUPRC values for each algorithm
            - AUROC: A dictionary containing AUROC values for each algorithm
    '''

    if userReferenceNetworkFile is None:
    # Read file for trueEdges
        trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                    '/' +dataDict['trueEdges'],
                                    sep = ',',
                                    header = 0, index_col = None)
    else:
        trueEdgesDF = pd.read_csv(str(userReferenceNetworkFile), sep = ',', header = 0, index_col = None)


    if onlyEdgesFromTFs:
        if tfsFile is None:
            print('ERROR: Please specify the path to a file containing a list of transcription factors in order to ignore edges from the transcription factors.')
        else:
            tfsDF = pd.read_csv(str(tfsFile), index_col = None, names=['Gene1'])
            trueEdgesDF = pd.merge(trueEdgesDF, tfsDF, on='Gene1', how='inner')

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

def computeScores(trueEdgesDF, predEdgeDF,
                  directed = True, selfEdges = True):
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

        selfEdges: bool
            A flag to indicate whether to include
            self-edges (selfEdges = True) or
            exclude self-edges (selfEdges = False) from evaluation.

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
        if selfEdges:
            possibleEdges = list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                         repeat = 2))
        else:
            possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                         r = 2))

        TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}

        # Compute TrueEdgeDict Dictionary
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth
        for key in TrueEdgeDict.keys():
            if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                   (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
                    TrueEdgeDict[key] = 1

        for key in PredEdgeDict.keys():
            subDF = predEdgeDF.loc[(predEdgeDF['Gene1'] == key.split('|')[0]) &
                               (predEdgeDF['Gene2'] == key.split('|')[1])]
            if len(subDF)>0:
                PredEdgeDict[key] = np.abs(subDF.EdgeWeight.values[0])

    # if undirected
    else:

        # Initialize dictionaries with all
        # possible edges
        if selfEdges:
            possibleEdges = list(combinations_with_replacement(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                                               r = 2))
        else:
            possibleEdges = list(combinations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                                               r = 2))
        TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}

        # Compute TrueEdgeDict Dictionary
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth

        for key in TrueEdgeDict.keys():
            if len(trueEdgesDF.loc[((trueEdgesDF['Gene1'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene2'] == key.split('|')[1])) |
                              ((trueEdgesDF['Gene2'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene1'] == key.split('|')[1]))]) > 0:
                TrueEdgeDict[key] = 1

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

    fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
                                     y_score=outDF['PredEdges'], pos_label=1)

    prec, recall, thresholds = precision_recall_curve(y_true=outDF['TrueEdges'],
                                                      probas_pred=outDF['PredEdges'], pos_label=1)

    return prec, recall, fpr, tpr, auc(recall, prec), auc(fpr, tpr)
