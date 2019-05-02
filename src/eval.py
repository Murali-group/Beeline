import pandas as pd
import numpy as np
from pathlib import Path
# for some reason these packages don't work on baobab
#import matplotlib.pyplot as plt
#import seaborn as sns
#sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from sklearn.metrics import precision_recall_curve, roc_curve, auc, roc_auc_score
from itertools import product, combinations_with_replacement
from tqdm import tqdm


def Eval(dataDict, inputSettings, runners):
    # Read file for trueEdges
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', 
                                header = 0, index_col = None)
    TrueEdgeDict = {'|'.join(p):0 for p in list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),repeat =2))}
    for key in TrueEdgeDict.keys():
        if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
            TrueEdgeDict[key] = 1
    # part 2 - Compute PR and ROC curves
    # by treating edges in the reference network as undirected
    uTrueEdgeDict = {'|'.join(p):0 for p in list(combinations_with_replacement(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),r = 2))}
    for key in uTrueEdgeDict.keys():
        if len(trueEdgesDF.loc[((trueEdgesDF['Gene1'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene2'] == key.split('|')[1])) |
                              ((trueEdgesDF['Gene2'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene1'] == key.split('|')[1]))])>0:
            uTrueEdgeDict[key] = 1

    # Initialize data dictionaries
    AUPRC, AUROC, FMAX = {}, {}, {}
    uAUPRC, uAUROC, uFMAX = {}, {}, {}
    num_edges = {}
    alg_name = {}
    for i, runner in tqdm(runners.items()):

        # check if the output rankedEdges file exists
        if Path(runner.final_ranked_edges).exists():
            predDF = pd.read_csv(
                runner.final_ranked_edges, sep='\t', header= 0, index_col=None)
            num_edges[runner.params_str] = len(predDF)
            alg_name[runner.params_str] = runner.name
            # if there are no edges, then give them all a value of 0
            if num_edges == 0:
                for d in [AUPRC, AUROC, FMAX, uAUPRC, uAUROC, uFMAX]:
                    d[runner.params_str] = 0
                continue

            auprc, auroc, fmax = compute_eval_measures(trueEdgesDF, TrueEdgeDict, predDF, directed=True)
            AUPRC[runner.params_str] = auprc
            AUROC[runner.params_str] = auroc
            FMAX[runner.params_str] = fmax

            auprc, auroc, fmax = compute_eval_measures(trueEdgesDF, uTrueEdgeDict, predDF, directed=False)
            uAUPRC[runner.params_str] = auprc
            uAUROC[runner.params_str] = auroc
            uFMAX[runner.params_str] = fmax
        else:
            print(runner.final_ranked_edges, ' does not exist. Skipping...')

    results = {'alg': alg_name, 'num_edges': num_edges, 'AUPRC': AUPRC, 'AUROC': AUROC, 'FMAX': FMAX,
                'uAUPRC': uAUPRC, 'uAUROC': uAUROC, 'uFMAX': uFMAX}
    return results


def compute_eval_measures(trueEdgesDF, TrueEdgeDict, predDF, directed=True):
    if directed:
        PredEdgeDict = {'|'.join(p):0 for p in list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),repeat =2))}
        for key in PredEdgeDict.keys():
            subDF = predDF.loc[(predDF['Gene1'] == key.split('|')[0]) &
                            (predDF['Gene2'] == key.split('|')[1])]
            if len(subDF)>0:
                PredEdgeDict[key] = np.abs(subDF.EdgeWeight.values[0])
    else:
        PredEdgeDict = {'|'.join(p):0 for p in list(combinations_with_replacement(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),r =2))}
        for key in PredEdgeDict.keys():
            subDF = predDF.loc[((predDF['Gene1'] == key.split('|')[0]) &
                                (predDF['Gene2'] == key.split('|')[1])) |
                                ((predDF['Gene2'] == key.split('|')[0]) &
                                (predDF['Gene1'] == key.split('|')[1]))]
            if len(subDF)>0:
                PredEdgeDict[key] = max(np.abs(subDF.EdgeWeight.values))

    outDF = pd.DataFrame([TrueEdgeDict,PredEdgeDict]).T
    outDF.columns = ['TrueEdges','PredEdges']
    #fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
    #                                    y_score=outDF['PredEdges'], pos_label=1)
    auroc = roc_auc_score(y_true=outDF['TrueEdges'], y_score=outDF['PredEdges'])

    prec, recall, thresholds = precision_recall_curve(y_true=outDF['TrueEdges'],
                                                        probas_pred=outDF['PredEdges'], pos_label=1)
    auprc = auc(recall, prec)
    #auroc = auc(fpr, tpr)
    fmax = compute_fmax(prec, recall)
    return auprc, auroc, fmax


def compute_fmax(prec, rec):
    f_measures = []
    for i in range(len(prec)):
        p, r = prec[i], rec[i]
        if p+r == 0:
            harmonic_mean = 0
        else:
            harmonic_mean = (2*p*r)/(p+r)
        f_measures.append(harmonic_mean)
    return max(f_measures)
