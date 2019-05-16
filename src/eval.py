import sys
import pandas as pd
import numpy as np
from pathlib import Path
# for some reason these packages don't work on baobab
#import matplotlib.pyplot as plt
#import seaborn as sns
#sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from sklearn.metrics import precision_recall_curve, roc_curve, auc, roc_auc_score
from itertools import product, combinations_with_replacement, permutations, combinations
from tqdm import tqdm


def Eval(dataDict, inputSettings, runners, plot_curves=False):
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    # Read file for trueEdges
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', 
                                header = 0, index_col = None)
    TrueEdgeDict = {'|'.join(p):0 for p in list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),r =2))}
    for key in TrueEdgeDict.keys():
        if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
            TrueEdgeDict[key] = 1
    # part 2 - Compute PR and ROC curves
    # by treating edges in the reference network as undirected
    uTrueEdgeDict = {'|'.join(p):0 for p in list(combinations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),r = 2))}
    for key in uTrueEdgeDict.keys():
        if len(trueEdgesDF.loc[((trueEdgesDF['Gene1'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene2'] == key.split('|')[1])) |
                              ((trueEdgesDF['Gene2'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene1'] == key.split('|')[1]))])>0:
            uTrueEdgeDict[key] = 1

    # Initialize data dictionaries
    AUPRC, AUROC, FMAX = {}, {}, {}
    uAUPRC, uAUROC, uFMAX = {}, {}, {}
    # for now, just do undirected
    #ePrec, r0_1Prec, r0_15Prec, r0_2Prec = {}, {}, {}, {}
    uePrec, ueRec = {}, {}
    ur0_1Prec, ur0_15Prec, ur0_2Prec = {}, {}, {} 
    num_edges = {}
    alg_name = {}
    recallDict, precisionDict, FPRDict, TPRDict = {}, {}, {}, {}
    for i, runner in tqdm(runners.items()):

        # check if the output rankedEdges file exists
        if Path(runner.final_ranked_edges).exists():
            tqdm.write("\tprocessing %s" % (runner.final_ranked_edges))
            predDF = pd.read_csv(
                runner.final_ranked_edges, sep='\t', header= 0, index_col=None)
            num_edges[runner.params_str] = len(predDF)
            alg_name[runner.params_str] = runner.name
            # if there are no edges, then give them all a value of 0
            if num_edges == 0:
                for d in [AUPRC, AUROC, FMAX, uAUPRC, uAUROC, uFMAX]:
                    d[runner.params_str] = 0
                continue
            key = runner.params_str

            if plot_curves is True:
                print("plot curves")
                key = runner.name
                # TODO, just use one parameter set for now
                tpr, fpr, prec, recall = compute_eval_measures(trueEdgesDF, TrueEdgeDict, predDF, directed=True, early=False, curves=True)
                recallDict[key] = recall
                precisionDict[key] = prec
                FPRDict[key] = fpr
                TPRDict[key] = tpr

            auprc, auroc, fmax = compute_eval_measures(trueEdgesDF, TrueEdgeDict, predDF, directed=True, early=False)
            AUPRC[key] = auprc
            AUROC[key] = auroc
            FMAX[key] = fmax
            #eauprc = compute_eval_measures(trueEdgesDF, TrueEdgeDict, predDF, directed=True, early=True)
            #eAUPRC[runner.params_str] = eauprc

            eprec, erec = compute_eval_measures(
                trueEdgesDF, uTrueEdgeDict, predDF, directed=False, early=True)
            uePrec[key], ueRec[key] = eprec, erec
            #eprec, (r0_1prec, r0_15prec, r0_2prec) = compute_eval_measures(
            #    trueEdgesDF, TrueEdgeDict, predDF, directed=True, early=True, recall_to_test=[0.1, 0.15, 0.2])
            #ePrec[key], r0_1Prec[key], r0_15Prec[key], r0_2Prec[key] = eprec, r0_1prec, r0_15prec, r0_2prec

            auprc, auroc, fmax = compute_eval_measures(trueEdgesDF, uTrueEdgeDict, predDF, directed=False)
            uAUPRC[key] = auprc
            uAUROC[key] = auroc
            uFMAX[key] = fmax
        else:
            print(runner.final_ranked_edges, ' does not exist. Skipping...')

    if plot_curves is True:
        plot_alg_curves(outDir, recallDict, precisionDict, FPRDict, TPRDict, AUPRC, AUROC)

    #results = {'alg': alg_name, 'num_edges': num_edges, 'AUPRC': AUPRC, 'eAUPRC': eAUPRC, 'AUROC': AUROC, 'FMAX': FMAX,
    results = {'alg': alg_name, 'num_edges': num_edges, 'AUPRC': AUPRC, 'AUROC': AUROC, 'FMAX': FMAX,
    #results = {'alg': alg_name, 'num_edges': num_edges, 
    'uePrec': uePrec, 'ueRec': ueRec,
                #'ePrec': ePrec, 'r0.1Prec': r0_1Prec, 'r0.15Prec': r0_15Prec, 'r0.2Prec': r0_2Prec}
    #results = {'alg': alg_name, 'num_edges': num_edges, 'AUPRC': AUPRC, 'AUROC': AUROC, 'FMAX': FMAX,
                'uAUPRC': uAUPRC, 'uAUROC': uAUROC, 'uFMAX': uFMAX
                }
    return results


def plot_alg_curves(outDir, recallDict, precisionDict, FPRDict, TPRDict, AUPRC, AUROC, directed=True):
    if directed:
        PRName = "/PRplot"
        ROCName = "/ROCplot"
    else:
        PRName = "/uPRplot"
        ROCName = "/uROCplot"
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


def get_subnetwork(trueEdgesDF, predDF, directed=True):
    # first limit the edges in the predDF to those present in the trueEdgesDF
    nodes = np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']])
    predDF = predDF[predDF['Gene1'].isin(nodes) & predDF['Gene2'].isin(nodes)]
    # remove self edges
    predDF = predDF[predDF['Gene1'] != predDF['Gene2']]
    # take the absolute value of the prediction scores
    predDF['EdgeWeight'] = predDF['EdgeWeight'].apply(np.abs)
    predDF = predDF.sort_values('EdgeWeight', ascending=False)
    predDF = predDF.reset_index()
    #print(len(predDF))
    # limit the predDF to the same size as the trueEdgesDF
    # keep ties
    num_edges = len(trueEdgesDF)
    if directed is False:
        num_edges = len(trueEdgesDF)*2
    #print(num_edges)
    score_at_num_edges = predDF.loc[num_edges-1]['EdgeWeight']
    predDF = predDF[predDF['EdgeWeight'] >= score_at_num_edges]
    # also, only get edges with a score > 0
    predDF = predDF[predDF['EdgeWeight'] > 0]

    return predDF


def compute_eval_measures(trueEdgesDF, TrueEdgeDict, predDF, directed=True, early=False, curves=False, recall_to_test=None):
    # TODO implement early undirected
    if early:
        predDF = get_subnetwork(trueEdgesDF, predDF, directed=directed)
    num_recovered = 0
    #print("True:")
    #print(len(trueEdgesDF))
    #print(trueEdgesDF)
    if directed:
        PredEdgeDict = {'|'.join(p):0 for p in list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),r =2))}
        for key in PredEdgeDict.keys():
            subDF = predDF.loc[(predDF['Gene1'] == key.split('|')[0]) &
                            (predDF['Gene2'] == key.split('|')[1])]
            if len(subDF)>0:
                PredEdgeDict[key] = np.abs(subDF.EdgeWeight.values[0])
                num_recovered += 1 
    else:
        PredEdgeDict = {'|'.join(p):0 for p in list(combinations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),r =2))}
        for key in PredEdgeDict.keys():
            subDF = predDF.loc[((predDF['Gene1'] == key.split('|')[0]) &
                                (predDF['Gene2'] == key.split('|')[1])) |
                                ((predDF['Gene2'] == key.split('|')[0]) &
                                (predDF['Gene1'] == key.split('|')[1]))]
            if len(subDF)>0:
                PredEdgeDict[key] = max(np.abs(subDF.EdgeWeight.values))

    #print(len(PredEdgeDict))
    outDF = pd.DataFrame([TrueEdgeDict,PredEdgeDict]).T
    outDF.columns = ['TrueEdges','PredEdges']
    #fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
    #                                    y_score=outDF['PredEdges'], pos_label=1)
    if early:
        # compute prec and rec, cutting off the curve at the recall of the # of edges
        prec, recall, thresholds = precision_recall_curve(
            y_true=outDF['TrueEdges'], probas_pred=outDF['PredEdges'], pos_label=1)
        # the second to last values are the early final prec/rec values 
        prec = prec[::-1]
        recall = recall[::-1]
        p, r = prec[-2], recall[-2]
        #auprc = auc(rec, prec)
        #auroc = auc(fpr, tpr)
        #fmax = compute_fmax(prec, rec)
        return p, r
    elif recall_to_test is not None:
        pass
    elif curves:
        fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'], y_score=outDF['PredEdges'])
        prec, recall, thresholds = precision_recall_curve(
            y_true=outDF['TrueEdges'], probas_pred=outDF['PredEdges'], pos_label=1)
        return tpr, fpr, prec, recall
    else:
        auroc = roc_auc_score(y_true=outDF['TrueEdges'], y_score=outDF['PredEdges'])

        prec, recall, thresholds = precision_recall_curve(
            y_true=outDF['TrueEdges'], probas_pred=outDF['PredEdges'], pos_label=1)
        auprc = auc(recall, prec)
        #auroc = auc(fpr, tpr)
        fmax = compute_fmax(prec, recall)
        return auprc, auroc, fmax


def compute_prec_rec(outDF):
    precision = [1]
    recall = [0]
    TP = 0
    FP = 0
    num_pos = len(outDF[outDF['TrueEdges'] == 1])
    outDF = outDF.sort_values('PredEdges', ascending=False)
    # get only the edges with a score > 0.
    outDF = outDF[outDF['PredEdges'] > 0]

    for label in outDF['TrueEdges'].values:
        if label == 1:
            TP += 1
            # precisions is the # of true positives / # true positives + # of false positives (or the total # of predictions)
            precision.append(TP / float(TP + FP))
            # recall is the # of recovered positives TP / TP + FN (total # of positives)
            rec = TP / float(num_pos)
            recall.append(rec)
        else:
            FP += 1
    return precision, recall

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
