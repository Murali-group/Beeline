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
        
        
    
def run(RunnerObj):
    '''
    Function to run GRISLI algorithm
    '''
    
    L = str(RunnerObj.params['L'])
    R = str(RunnerObj.params['R'])
    alphaMin = str(RunnerObj.params['alphaMin'])
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/GRISLI/"
    os.makedirs(outDir, exist_ok = True)
    
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        inputPath = "data"+str(RunnerObj.inputDir).split(str(Path.cwd()))[1]+"/GRISLI/"+str(idx)+"/"
        os.makedirs(outDir+str(idx), exist_ok = True)

        outFile = "data/" +  str(outDir) +str(idx)+"/outFile.txt"

        cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/runGRISLI/data/ grisli:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time'+str(idx)+'.txt', './GRISLI ',inputPath, outFile, L, R, alphaMin,'\"'])
    
        print(cmdToRun)
        os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from GRISLI.
    '''
    
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/GRISLI/"

    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
    colNames = PTData.columns
    OutSubDF = [0]*len(colNames)
    
    for indx in range(len(colNames)):
        # Read output
        outFile = str(indx)+'/outFile.txt'
        if not Path(outDir+outFile).exists():
            # Quit if output file does not exist
            print(outDir+outFile+' does not exist, skipping...')
            return
        OutDF = pd.read_csv(outDir+outFile, sep = ',', header = None)    
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
        outFileName = outDir + str(indx)+ '/rankedEdges.csv'
        outFile = open(outFileName,'w')
        outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

        for row, col, val in zip(rows, cols, DFSorted):
            outFile.write('\t'.join([GeneList[row],GeneList[col],str((len(GeneList)*len(GeneList))-val)])+'\n')
        outFile.close()
        
        OutSubDF[indx] = pd.read_csv(outFileName, sep = '\t', header = 0)

        # megre the dataframe by taking the maximum value from each DF
        # From here: https://stackoverflow.com/questions/20383647/pandas-selecting-by-label-sometimes-return-series-sometimes-returns-dataframe
    outDF = pd.concat(OutSubDF)

    res = outDF.groupby(['Gene1','Gene2'],as_index=False).max()
    #print(res.head())
    # Sort values in the dataframe   
    finalDF = res.sort_values('EdgeWeight',ascending=False)  
    
    finalDF.to_csv(outDir+'rankedEdges.csv',sep='\t', index = False)
