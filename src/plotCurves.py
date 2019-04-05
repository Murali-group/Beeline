import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")

def EvalCurves(dataDict, inputSettings):
    '''
    Computes PR and ROC curves
    for a given dataset and a set of algorithms
    '''
    
    # Read file for trueEdges
    trueEdgesFile = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = '\t', header = 0, index_col = None)


    
    # Initialize data dictionaries
    precisionDict = {}
    recallDict = {}
    FPRDict = {}
    AUPRC = {}
    AUROC = {}
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    for algo in inputSettings.algorithms:

        # check if the output rankedEdges file exists
        if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
            precisionDict[algo[0]] = [] # Initialize Precsion
            recallDict[algo[0]] = [] # Initialize Recall
            FPRDict[algo[0]] = [] # Initialize FPR
            predEdgesFile = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                        sep = '\t', header =  0, index_col = None)
            
            
            tp = 0
            fp = 0
            total = 0 # total predictions made
            totalTrue = trueEdgesFile.shape[0] # Condition Positives
            pOld = 0
            rOld = 0
            AUPRC[algo[0]] = 0 # Initialize AUPRC
            AUROC[algo[0]] = 0 # Initialize AUROC
            
            for idx, row in predEdgesFile.iterrows():
                if trueEdgesFile.loc[(trueEdgesFile['Gene1'] == row['Gene1']) & \
                                     (trueEdgesFile['Gene2'] == row['Gene2'])].shape[0] > 0:
                    tp += 1
                else:
                    fp += 1
                total += 1
                
                pNew = float(tp)/float(total)
                rNew = float(tp)/float(totalTrue)
                
                precisionDict[algo[0]].append(pNew)
                recallDict[algo[0]].append(rNew)
                
                AUPRC[algo[0]] += ((rNew - rOld)*(pOld + pNew)/2) # compute AUPRC
                
                pOld = pNew
                rOld = rNew
                FPRDict[algo[0]].append(float(fp)) # List of FP values
                
            FPRDict[algo[0]] = [val/float(total - totalTrue) for val in FPRDict[algo[0]]] # update FPR
            tprOld = 0
            fprOld = 0
            for idx in range(len(FPRDict[algo[0]])):
                tprNew = recallDict[algo[0]][idx]
                fprNew = FPRDict[algo[0]][idx]
                AUROC[algo[0]] += ((fprNew - fprOld)*(tprOld + tprNew)/2)
                tprOld = tprNew
                fprOld = fprNew
        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' does not exist. Skipping...')
            
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
        sns.lineplot(FPRDict[key],recallDict[key], ci=None)
        legendList.append(key + ' (AUROC = ' + str("%.2f" % (AUROC[key]))+')')
        
    plt.plot([0, 1], [0, 1], linewidth = 1.5, color = 'k', linestyle = '--')

    plt.xlim(0,1)    
    plt.ylim(0,1)
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.legend(legendList) 
    plt.savefig(outDir+'/ROCplot.pdf')
    plt.savefig(outDir+'/ROCplot.png')
