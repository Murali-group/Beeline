import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for TENET.
    If the folder/files under RunnerObj.datadir exist,
    this function will not do anything.
    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    if not RunnerObj.inputDir.joinpath("TENET").exists():
        print("Input folder for TENET does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("TENET").mkdir(exist_ok = False)
            
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        RunnerObj.inputDir.joinpath("TENET/"+str(idx)).mkdir(exist_ok = True)
        
        # Select cells belonging to each pseudotime trajectory
        colName = colNames[idx]
        index = PTData[colName].index[PTData[colName].notnull()]
        
        exprName = "TENET/"+str(idx)+"/ExpressionData.csv"
        ExpressionData.loc[:,index].T.to_csv(RunnerObj.inputDir.joinpath(exprName),
                                 sep = ',', header  = False, index = False)
        
        cellName = "TENET/"+str(idx)+"/PseudoTime.csv"
        ptDF = PTData.loc[index,[colName]]                
        ptDF.to_csv(RunnerObj.inputDir.joinpath(cellName),
                                 sep = ',', header  = False, index = False)
    
    # TODO: Create the cell-selection raw text file for the TENET algorithm
    # if not RunnerObj.inputDir.joinpath("TENET/cell_select.txt"):
    #     pass

def run(RunnerObj):
    
    
    # # tenet wants info as gene columns and cell rows, opposite of how info is provided through pipeline
    # toTranspose = pd.read_csv(inputPath)
    # transposedDF = pd.T
    # #remove nontransposed version
    # os.remove(inputPath)
    # transposedDF.to_csv(inputPath)
    
    
    # make output dirs if they do not exist:
    
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/TENET/"
    cellHistoryLength = str(RunnerObj.params['historyLength'])
    threadsToRun = str(RunnerObj.params['threads'])
    PTData = pd.read_csv(RunnerObj.inputDir.joinPath(RunnerObj.cellData), header=0, index_col=0)
    colNames = PTData.columns
    outPath = "data/" +  str(outDir) + 'outFile.txt'

    for idx in range(len(colNames)):
        expressionPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/TENET/" + str(idx) + "/ExpressionData.csv" 
        PTPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/TENET/" + str(idx) + "/PsuedoTime.csv"
        os.makedirs(outDir, exist_ok = True)
        cellNum = pd.read_csv(expressionPath).shape[0]
        cmdToRun = ' '.join(['docker run --rm -v', 
                            str(Path.cwd())+':/data tenet:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time.txt',
                            'yes 1 | head -n ' + cellNum + ' > cellSelect.txt;',
                            '/TENET',
                            expressionPath,
                            threadsToRun,
                            PTPath,
                            'cellSelect.txt',
                            cellHistoryLength])
    

def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCSGL
    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    pass