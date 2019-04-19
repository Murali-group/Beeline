import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for SCRIBE.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("SCRIBE").exists():
        print("Input folder for SCRIBE does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("SCRIBE").mkdir(exist_ok = False)
    
    
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath('ExpressionData.csv'),
                                     header = 0, index_col = 0)
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath('PseudoTime.csv'),
                             header = 0, index_col = 0)
    
    if not RunnerObj.inputDir.joinpath("SCRIBE/ExpressionData.csv").exists(): 
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath("SCRIBE/ExpressionData.csv"),
                             sep = ',', header  = True, index = True)
        
    if not RunnerObj.inputDir.joinpath("SCRIBE/CellData.csv").exists():
        del PTData.index.name
        PTData.to_csv(RunnerObj.inputDir.joinpath("SCRIBE/CellData.csv"),
                             sep = ',', header  = True, index = True)
        
    if not RunnerObj.inputDir.joinpath("SCRIBE/GeneData.csv").exists():
        # required column!!
        geneDict = {}
        geneDict['gene_short_name'] = [gene.replace('x_', '') for gene in ExpressionData.index]
        
        geneDF = pd.DataFrame(geneDict, index = ExpressionData.index)
        geneDF.to_csv(RunnerObj.inputDir.joinpath("SCRIBE/GeneData.csv"), 
                      sep = ',', header = True)
    
def run(RunnerObj):
    '''
    Function to run SCRIBE algorithm.
    To see all the inputs runScribe.R script takes, run:
    docker run scribe:base /bin/sh -c "Rscript runScribe.R -h"
    '''
    
    inputPath = "data/"+str(RunnerObj.inputDir).split("ModelEval/")[1]+"/SCRIBE/"

    
    # required inputs
    delay = str(RunnerObj.params['delay'])
    method = str(RunnerObj.params['method'])
    low = str(RunnerObj.params['lowerDetectionLimit'])
    fam = str(RunnerObj.params['expressionFamily'])

    # optional inputs
    log = str(RunnerObj.params['log'])
    ignorePT = str(RunnerObj.params['ignorePT'])
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SCRIBE/"
    os.makedirs(outDir, exist_ok = True)

    # Build the command to run Scribe
    
    cmdToRun = ' '.join(['docker run --rm -v ~/ModelEval:/data/ scribe:base /bin/sh -c \"Rscript runScribe.R',
                   '-e',inputPath +'ExpressionData.csv', '-c',inputPath + 'CellData.csv', 
                   '-g',inputPath + 'GeneData.csv', '-o data/'+outDir, '-d',delay, '-l', low,
                   '-m', method, '-x',fam])
    
    if RunnerObj.params['log'] == 'True':
        cmdToRun += ' --log'
    if RunnerObj.params['ignorePT'] == 'True':
        cmdToRun += ' --i'
        
    cmdToRun += '\"'
    
    print(cmdToRun)
    sys.exit()
    os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCRIBE.
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SCRIBE/"
    if not Path(outDir).exists():
        raise FileNotFoundError()
        
    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = ',', header = 0, index_col = 0)
    
    # Sort values in a matrix using code from:
    # https://stackoverflow.com/questions/21922806/sort-values-of-matrix-in-python
    
    OutDF = OutDF.clip(lower=0) # if less than 0, set it to zero
    print(OutDF.head())
    idx = np.argsort(OutDF.values, axis = None)[::-1]
    print(idx)
    rows, cols = np.unravel_index(idx, OutDF.shape)    
    DFSorted = OutDF[rows, cols]
    
    # read input file for list of gene names
    GeneList = list(OutDF.index)
    
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for row, col, val in zip(rows, cols, DFSorted):
        outFile.write('\t'.join([GeneList[row],GeneList[col],str(val)])+'\n')
    outFile.close()
    