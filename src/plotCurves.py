import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from itertools import product, combinations_with_replacement


def EvalCurves(dataDict, inputSettings):
    '''
    Computes PR and ROC curves
    for a given dataset and a set of algorithms
    '''
    
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
    ensembleFinal.sort_values(by='EdgeWeight',ascending = True, inplace = True)
    print("Writing Ensemble Model output to ",outDir)
    ensembleFinal[['Gene1','Gene2','EdgeWeight']].to_csv(outDir+'/rankedEdges.csv', sep = '\t', index = False)            


            
    # Initialize data dictionaries
    precisionDict = {}
    recallDict = {}
    FPRDict = {}
    TPRDict = {}
    AUPRC = {}
    AUROC = {}
    for algo in inputSettings.algorithms:

        # check if the output rankedEdges file exists
        if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
             # Initialize Precsion

            
            predDF = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                        sep = '\t', header =  0, index_col = None)

            PredEdgeDict = {'|'.join(p):0 for p in list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),repeat =2))}
            for key in PredEdgeDict.keys():
                subDF = predDF.loc[(predDF['Gene1'] == key.split('|')[0]) &
                                   (predDF['Gene2'] == key.split('|')[1])]
                if len(subDF)>0:
                    PredEdgeDict[key] = np.abs(subDF.EdgeWeight.values[0])
            
            outDF = pd.DataFrame([TrueEdgeDict,PredEdgeDict]).T
            outDF.columns = ['TrueEdges','PredEdges']
            fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
                                             y_score=outDF['PredEdges'], pos_label=1)
            
            prec, recall, thresholds = precision_recall_curve(y_true=outDF['TrueEdges'],
                                                              probas_pred=outDF['PredEdges'], pos_label=1)
            precisionDict[algo[0]] = prec
            recallDict[algo[0]] = recall 
            FPRDict[algo[0]] = fpr 
            TPRDict[algo[0]] = tpr
            AUPRC[algo[0]] = auc(recall, prec)
            AUROC[algo[0]] = auc(fpr, tpr)
        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' does not exist. Skipping...')
            
            
    predDF = pd.read_csv(outDir + '/rankedEdges.csv', \
                                sep = '\t', header =  0, index_col = None)

    PredEdgeDict = {'|'.join(p):0 for p in list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),repeat =2))}
    for key in PredEdgeDict.keys():
        subDF = predDF.loc[(predDF['Gene1'] == key.split('|')[0]) &
                           (predDF['Gene2'] == key.split('|')[1])]
        if len(subDF)>0:
            PredEdgeDict[key] = subDF.EdgeWeight.values[0]

    outDF = pd.DataFrame([TrueEdgeDict,PredEdgeDict]).T
    outDF.columns = ['TrueEdges','PredEdges']
    fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
                                     y_score=outDF['PredEdges'], pos_label=1)

    prec, recall, thresholds = precision_recall_curve(y_true=outDF['TrueEdges'],
                                                      probas_pred=outDF['PredEdges'], pos_label=1)
    precisionDict['Ensemble'] = prec
    recallDict['Ensemble'] = recall 
    FPRDict['Ensemble'] = fpr 
    TPRDict['Ensemble'] = tpr
    AUPRC['Ensemble'] = auc(recall, prec)
    AUROC['Ensemble'] = auc(fpr, tpr)

    ## Make PR curves
    legendList = []
    for key in precisionDict.keys():
        print(key)
        sns.lineplot(recallDict[key],precisionDict[key], ci=None)
        legendList.append(str(key) + ' (AUPRC = ' + str("%.2f" % (AUPRC[key]))+')')

    plt.xlim(0,1)    
    plt.ylim(0,1)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(legendList) 
    plt.savefig(outDir+'/PRplot.pdf')
    plt.savefig(outDir+'/PRplot.png')
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
    plt.savefig(outDir+'/ROCplot.pdf')
    plt.savefig(outDir+'/ROCplot.png')
    plt.clf()

    # part 2 - Compute PR and ROC curves
    # by treating edges in the reference network
    # as undirected

    
    TrueEdgeDict = {'|'.join(p):0 for p in list(combinations_with_replacement(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),r = 2))}
    for key in TrueEdgeDict.keys():
        if len(trueEdgesDF.loc[((trueEdgesDF['Gene1'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene2'] == key.split('|')[1])) |
                              ((trueEdgesDF['Gene2'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene1'] == key.split('|')[1]))])>0:
            TrueEdgeDict[key] = 1
            
    # Initialize data dictionaries
    precisionDict = {}
    recallDict = {}
    FPRDict = {}
    TPRDict = {}
    uAUPRC = {}
    uAUROC = {}
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    for algo in inputSettings.algorithms:

        # check if the output rankedEdges file exists
        if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
             # Initialize Precsion

            predDF = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                        sep = '\t', header =  0, index_col = None)

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

            fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
                                             y_score=outDF['PredEdges'], pos_label=1)
            
            prec, recall, thresholds = precision_recall_curve(y_true=outDF['TrueEdges'],
                                                              probas_pred=outDF['PredEdges'], pos_label=1)
            precisionDict[algo[0]] = prec
            recallDict[algo[0]] = recall 
            FPRDict[algo[0]] = fpr 
            TPRDict[algo[0]] = tpr
            uAUPRC[algo[0]] = auc(recall, prec)
            uAUROC[algo[0]] = auc(fpr, tpr)
        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' does not exist. Skipping...')
            
    ## Make PR curves
    legendList = []
    for key in precisionDict.keys():
        print(key)
        sns.lineplot(recallDict[key],precisionDict[key], ci=None)
        legendList.append(str(key) + ' (AUPRC = ' + str("%.2f" % (uAUPRC[key]))+')')

    plt.xlim(0,1)    
    plt.ylim(0,1)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(legendList) 
    plt.savefig(outDir+'/uPRplot.pdf')
    plt.savefig(outDir+'/uPRplot.png')
    plt.clf()
    
    ## Make ROC curves
    legendList = []
    for key in recallDict.keys():
        sns.lineplot(FPRDict[key],TPRDict[key], ci=None)
        legendList.append(key + ' (AUROC = ' + str("%.2f" % (uAUROC[key]))+')')
        
    plt.plot([0, 1], [0, 1], linewidth = 1.5, color = 'k', linestyle = '--')

    plt.xlim(0,1)    
    plt.ylim(0,1)
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.legend(legendList) 
    plt.savefig(outDir+'/uROCplot.pdf')
    plt.savefig(outDir+'/uROCplot.png')
    plt.clf()
    return AUPRC, AUROC, uAUPRC, uAUROC
