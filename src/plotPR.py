import pandas as pd
import seaborn as sns
from pathlib import Path
import seaborn as sns

def PRCurves(dataDict, inputSettings):
    
    
    trueEdgesFile = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = '\t', header = 0, index_col = None)
    tp = 0
    fp = 0
    total = 0
    totalTrue = trueEdgesFile.shape[0]
    precision = []
    recall = []
    
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    for algo in inputSettings.algorithms:
        # check if the output rankedEdges file exists
        if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
            predEdgesFile = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                        sep = '\t', header =  0, index_col = None)
            
            for idx, row in predEdgesFile.iterrows():
                if trueEdgesFile.loc[(trueEdgesFile['Gene1'] == row['Gene1']) & \
                                     (trueEdgesFile['Gene2'] == row['Gene2'])].shape[0] > 0:
                    tp += 1
                else:
                    fp += 1
                total += 1
                precision.append(float(tp)/float(total))
                recall.append(float(tp)/float(totalTrue))
            sns_plot = sns.lineplot(recall, precision)
            sns_plot.figure.savefig('output.png')
        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' does not exist. Skipping...')
            
