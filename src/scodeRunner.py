import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for SCODE.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("SCODE").exists():
        print("Input folder for SCODE does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("SCODE").mkdir(exist_ok = False)
        

    
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        # Select cells belonging to each pseudotime trajectory
        colName = colNames[idx]
        index = PTData[colName].index[PTData[colName].notnull()]
        exprName = "SCODE/ExpressionData"+str(idx)+".csv"
        ExpressionData.loc[:,index].to_csv(RunnerObj.inputDir.joinpath(exprName),
                                 sep = '\t', header  = False, index = False)
        cellName = "SCODE/PseudoTime"+str(idx)+".csv"
        ptDF = PTData.loc[index,[colName]]        
        # SCODE expects a column labeled PseudoTime.
        ptDF.rename(columns = {colName:'PseudoTime'}, inplace = True)
        # output file
        ptDF.to_csv(RunnerObj.inputDir.joinpath(cellName),
                                 sep = '\t', header  = False)
        
    
def run(RunnerObj):
    '''
    Function to run SCODE algorithm
    '''
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SCODE/"
    os.makedirs(outDir, exist_ok = True)
    
    inputPath = "data"+str(RunnerObj.inputDir).split(str(Path.cwd()))[1]+"/SCODE/"

    z = str(RunnerObj.params['z'])
    nIter = str(RunnerObj.params['nIter'])
    nRep = str(RunnerObj.params['nRep'])
    
    
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    

    for idx in range(len(colNames)):

        ExpressionData = pd.read_csv(str(RunnerObj.inputDir).split(str(Path.cwd())+'/')[1]+"/SCODE/"+\
                                     "ExpressionData"+str(idx)+".csv",
                                     header = None, index_col = None, sep ='\t')
        nCells = str(ExpressionData.shape[1])
        nGenes = str(ExpressionData.shape[0])

        os.makedirs(outDir+str(idx), exist_ok = True)

        cmdToRun = ' '.join(['docker run --rm -v', 
                             str(Path.cwd())+':/SCODE/data/  scode:base /bin/sh -c \"time -v -o',
                             "data/" + str(outDir) + 'time'+str(idx)+'.txt', 'ruby run_R.rb',
                            inputPath +'ExpressionData'+str(idx)+'.csv', 
                            inputPath + 'PseudoTime'+str(idx)+'.csv', 
                            'data/'+outDir+str(idx),
                             nGenes, z, nCells, nIter, nRep, '\"'])
        print(cmdToRun)
        os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCODE.
    '''
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SCODE/"

    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
    colNames = PTData.columns
    for indx in range(len(colNames)):
        # Read output
        outFile = str(indx)+'/meanA.txt'
        if not Path(outDir+outFile).exists():
            # Quit if output file does not exist

            print(outDir+outFile+' does not exist, skipping...')
            return
        OutDF = pd.read_csv(outDir+outFile, sep = '\t', header = None)



        # Sort values in a matrix using code from:
        # https://stackoverflow.com/questions/21922806/sort-values-of-matrix-in-python
        OutMatrix = np.abs(OutDF.values)
        idx = np.argsort(OutMatrix, axis = None)[::-1]
        rows, cols = np.unravel_index(idx, OutDF.shape)    
        DFSorted = OutMatrix[rows, cols]

        # read input file for list of gene names
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                         header = 0, index_col = 0)
        GeneList = list(ExpressionData.index)

        outFile = open(outDir + 'outFile'+str(indx)+'.csv','w')
        outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

        for row, col, val in zip(rows, cols, DFSorted):
            outFile.write('\t'.join([GeneList[row],GeneList[col],str(val)])+'\n')
        outFile.close()
        
    OutSubDF = [0]*len(colNames)
    for indx in range(len(colNames)):
        outFile = 'outFile'+str(indx)+'.csv'
        OutSubDF[indx] = pd.read_csv(outDir+outFile, sep = '\t', header = 0)

        OutSubDF[indx].EdgeWeight = np.abs(OutSubDF[indx].EdgeWeight)

    outDF = pd.concat(OutSubDF)
    FinalDF = outDF[outDF['EdgeWeight'] == outDF.groupby(['Gene1','Gene2'])['EdgeWeight'].transform('max')]
    FinalDF.sort_values(['EdgeWeight'], ascending = False, inplace = True)
    FinalDF.to_csv(outDir+'rankedEdges.csv',sep = '\t', index = False)
