import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for GRISLI.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("GRISLI").exists():
        print("Input folder for GRISLI does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("GRISLI").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("GRISLI/ExpressionData.tsv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath('ExpressionData.csv'),
                                     header = 0, index_col = 0)
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath("GRISLI/ExpressionData.tsv"),
                             sep = '\t', header  = False, index = False)
        
    if not RunnerObj.inputDir.joinpath("GRISLI/PseudoTime.tsv").exists():
        PTData = pd.read_csv(RunnerObj.inputDir.joinpath('PseudoTime.csv'),
                             header = 0, index_col = 0)
        PTData['PseudoTime'].to_csv(RunnerObj.inputDir.joinpath("GRISLI/PseudoTime.tsv"),
                             sep = '\t', header  = False, index = False)
        
    
def run(RunnerObj):
    '''
    Function to run GRISLI algorithm
    '''
    inputPath = "data/"+str(RunnerObj.inputDir).split("RNMethods/")[1]+"/GRISLI/"
    
    L = str(RunnerObj.params['L'])
    R = str(RunnerObj.params['R'])
    alphaMin = str(RunnerObj.params['alphaMin'])
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/GRISLI/"
    os.makedirs(outDir, exist_ok = True)
    
    outFile = "data/" +  str(outDir) + 'outFile.txt'
    cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/runGRISLI/data/ grisli:test /bin/sh -c \"./GRISLI ', 
                         inputPath, outFile, L, R, alphaMin,'\"'])
    
    print(cmdToRun)
    os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from GRISLI.
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/GRISLI/"
    if not Path(outDir).exists():
        raise FileNotFoundError()
        
    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = ',', header = None)
    
    # Sort values in a matrix using code from:
    # https://stackoverflow.com/questions/21922806/sort-values-of-matrix-in-python
    OutMatrix = OutDF.values
    idx = np.argsort(OutMatrix, axis = None)
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
    