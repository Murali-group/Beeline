import os
import subprocess
import pandas as pd
from pathlib import Path
import numpy as np


params_order = ['L', 'R', 'alphaMin']
default_params = {'L': '10', 'R': '1500', 'alphaMin': '0.3'}

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for GRISLI.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("GRISLI").exists():
        print("Input folder for GRISLI does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("GRISLI").mkdir(exist_ok = False)
        
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                         header = 0, index_col = 0)
        PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                                 header = 0, index_col = 0)

        colNames = PTData.columns
        for idx in range(len(colNames)):
            RunnerObj.inputDir.joinpath("GRISLI/"+str(idx)).mkdir(exist_ok = True)
            
            # Select cells belonging to each pseudotime trajectory
            colName = colNames[idx]
            index = PTData[colName].index[PTData[colName].notnull()]
            
            exprName = "GRISLI/"+str(idx)+"/ExpressionData.tsv"
            ExpressionData.loc[:,index].to_csv(RunnerObj.inputDir.joinpath(exprName),
                                     sep = '\t', header  = False, index = False)
            
            cellName = "GRISLI/"+str(idx)+"/PseudoTime.tsv"
            ptDF = PTData.loc[index,[colName]]                
            ptDF.to_csv(RunnerObj.inputDir.joinpath(cellName),
                                     sep = '\t', header  = False, index = False)

    setupParams(RunnerObj)
        

def setupParams(RunnerObj):
    # if the parameters aren't specified, then use default parameters
    # TODO allow passing in multiple sets of hyperparameters
    params = RunnerObj.params
    for param, val in default_params.items():
        if param not in params:
            params[param] = val
    params_str = build_params_str(params)
    RunnerObj.params_str = params_str
    RunnerObj.params = params

    # the final file is written here:
    outDir = "outputs/" + str(RunnerObj.inputDir).split("inputs/")[1]+"/GRISLI/"
    print(outDir)
    RunnerObj.outDir = outDir
    RunnerObj.final_ranked_edges = "%s/%s-rankedEdges.csv" % (outDir, RunnerObj.params_str)


def build_params_str(params):
    params_str = "L%s-R%s-a%s" % (params['L'], params['R'], params['alphaMin'])
    return params_str


def run(RunnerObj):
    '''
    Function to run GRISLI algorithm
    '''
    inputPath = "data/"+str(RunnerObj.inputDir).split("RNMethods/")[1]+"/GRISLI/"

    L = str(RunnerObj.params['L'])
    R = str(RunnerObj.params['R'])
    alphaMin = str(RunnerObj.params['alphaMin'])

    # make output dirs if they do not exist:
    os.makedirs(RunnerObj.outDir, exist_ok = True)
    
    #RunnerObj.outFile = "data/" + str(RunnerObj.outDir) + RunnerObj.params_str + '-outFile.txt'
    RunnerObj.outFile = str(RunnerObj.outDir) + RunnerObj.params_str + '-outFile.txt'
    outFile = "data/" + RunnerObj.outFile
    print(RunnerObj.outFile)

    cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/runGRISLI/data/ grisli:base /bin/sh -c \"time -v -o', "data/" + str(RunnerObj.outDir) + 'time.txt', './GRISLI ',inputPath, outFile, L, R, alphaMin,'\"'])
    
    print(cmdToRun)
    #os.system(cmdToRun)
    subprocess.check_call(cmdToRun, shell=True)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from GRISLI.
    '''
    # Quit if output directory does not exist
    outFile = RunnerObj.outFile
    if not Path(outFile).exists():
        print(outFile+' does not exist, skipping...')
        return

    # Read output
    OutDF = pd.read_csv(outFile, sep = ',', header = None)

    # Sort values in a matrix using code from:
    # https://stackoverflow.com/questions/21922806/sort-values-of-matrix-in-python
    OutMatrix = OutDF.values
    idx = np.argsort(OutMatrix, axis = None)
    rows, cols = np.unravel_index(idx, OutDF.shape)    
    DFSorted = OutMatrix[rows, cols]

    # read input file for list of gene names
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    GeneList = list(ExpressionData.index)

    print("\twriting processed rankedEdges to %s" % (RunnerObj.final_ranked_edges))
    with open(RunnerObj.final_ranked_edges,'w') as out:
        out.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

        for row, col, val in zip(rows, cols, DFSorted):
            out.write('\t'.join([GeneList[row],GeneList[col],str((len(GeneList)*len(GeneList))-val)])+'\n')

