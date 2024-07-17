import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for DEEPSEM.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("DEEPSEM").exists():
        print("Input folder for DEEPSEM does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("DEEPSEM").mkdir(exist_ok = False)
        
    # TODO REMOVE COMMENT if not RunnerObj.inputDir.joinpath("DEEPSEM/ExpressionData.csv").exists():
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                    header = 0, index_col = 0)

    # Write .csv file
    ExpressionData.T.to_csv(RunnerObj.inputDir.joinpath("DEEPSEM/ExpressionData.csv"),
                             sep = ',', header  = True, index = True)
    # TODO REMOVE COMMENT if not RunnerObj.inputDir.joinpath("DEEPSEM/refNetwork.csv").exists():
    refNetwork = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges), 
                                header = 0, index_col = 0)
    refNetwork.to_csv(RunnerObj.inputDir.joinpath("DEEPSEM/refNetwork.csv"),
                            sep = ',', header  = True, index = True)
    
def run(RunnerObj):
    '''
    Function to run DEEPSEM algorithm
    '''
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/DEEPSEM/ExpressionData.csv"
    refPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/DEEPSEM/refNetwork.csv"
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/DEEPSEM/"
    os.makedirs(outDir, exist_ok = True)

    cellTypeSpecific = RunnerObj.params['isCellTypeSpecific']
    alpha = str(RunnerObj.params['alpha'])
    beta = str(RunnerObj.params['beta'])
    n_epochs = str(RunnerObj.params['numberEpochs'])

    if cellTypeSpecific == False:
        c_type = "--task=non_celltype_GRN"
    else:
        c_type = "--task=celltype_GRN"

    outPath = "data/" +  str(outDir) + 'outFile.tsv'
    cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/data/ --expose=41269', 
                         'deepsem:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time.txt', 
                         #'python -m pip list',
                         'python main.py',
                         c_type,
                         '--data_file='+inputPath,
                         '--net_file='+refPath,  
                         '--setting=new',
                         '--alpha='+alpha,
                         '--beta='+beta,
                         '--n_epoch='+n_epochs,
                         '--save_name='+outPath, 
                         '\"'])
    print(cmdToRun)
    os.system(cmdToRun)


def parseOutput(RunnerObj):
    pass
    '''
    Function to parse outputs from DEEPSEM.
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/DEEPSEM/"
    
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = '\t', header = 0)
    
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    # for idx, row in OutDF.iterrows():
    #     outFile.write('\t'.join([row['TF'],row['target'],str(row['importance'])])+'\n')
    outFile.close()
    
