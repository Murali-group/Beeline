import os
import pandas as pd
from pathlib import Path
import numpy as np
from sklearn import preprocessing

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for JUMP3.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("JUMP3").exists():
        print("Input folder for JUMP3 does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("JUMP3").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("JUMP3/ExpressionData.csv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        newExpressionData = ExpressionData.T.copy()
        PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
        # make sure the indices are strings for both dataframes
        newExpressionData.index = newExpressionData.index.map(str) 
        PTData.index = PTData.index.map(str) 
        # Acc. to JUMP3:
        # In input argument Time, the first time point of each time series must be 0.
        # Also has to be an integer!
        newExpressionData['Time'] = PTData['PseudoTime']-PTData['PseudoTime'].min()
        if 'Experiment' in PTData:
            newExpressionData['Experiment'] = PTData['Experiment']
        else:
            # generate it from cell number Ex_y, where x is experiment number
            #newExpressionData['Experiment'] = [int(x.split('_')[0].strip('E')) for x in PTData.index.astype(str)]
            newExpressionData['Experiment'] = 1
            
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath("JUMP3/ExpressionData.csv"),
                             sep = ',', header  = True, index = False)
    
    
def run(RunnerObj):
    '''
    Function to run GRN-VBEM algorithm
    '''
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/JUMP3/ExpressionData.csv"
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/JUMP3/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath = "data/" +  str(outDir) + 'outFile.txt'
    cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/JUMP3/data/ jump3:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time.txt', './runJump3',
                         inputPath, outPath, '\"'])
    print(cmdToRun)
    os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from JUMP3.
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/JUMP3/"
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
        
    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = ',')
    
    # Sort values in a matrix using code from:
    # https://stackoverflow.com/questions/21922806/sort-values-of-matrix-in-python
    OutMatrix = np.abs(OutDF.values)
    idx = np.argsort(OutMatrix, axis = None)[::-1]
    rows, cols = np.unravel_index(idx, OutDF.shape)    
    DFSorted = OutMatrix[rows, cols]
    
    # read input file for list of gene names
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath('ExpressionData.csv'),
                                     header = 0, index_col = 0)
    GeneList = list(ExpressionData.index)
    
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for row, col, val in zip(rows, cols, DFSorted):
        outFile.write('\t'.join([GeneList[row],GeneList[col],str(val)])+'\n')
    outFile.close()
    
