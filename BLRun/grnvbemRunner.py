import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for GRNVBEM.
    It will create the input folder at RunnerObj.dataset.name/GRNVBEM/ if it
    does not exist already. The input folder will contain an ExpressionData.csv with
    cells ordered according to the pseudotime along the columns, and genes along
    the rows. If the files already exist, this function will overwrite it.
    
    '''
    if not RunnerObj.inputDir.joinpath("GRNVBEM").exists():
        print("Input folder for GRNVBEM does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("GRNVBEM").mkdir(exist_ok = False)
        
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        # Select cells belonging to each pseudotime trajectory
        colName = colNames[idx]
        index = PTData[colName].index[PTData[colName].notnull()]
        exprName = "GRNVBEM/ExpressionData"+str(idx)+".csv"
        
        subPT = PTData.loc[index,:]
        subExpr = ExpressionData[index]
        # Order columns by PseudoTime
        newExpressionData = subExpr[subPT.sort_values([colName]).index.astype(str)]
        
        newExpressionData.insert(loc = 0, column = 'GENES', \
                                                     value = newExpressionData.index)


        # Write .csv file
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath(exprName),
                             sep = ',', header  = True, index = False)
    
def run(RunnerObj):
    '''
    Function to run GRN-VBEM algorithm
    '''
    
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1]
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/GRNVBEM/"
    os.makedirs(outDir, exist_ok = True)
    
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        exprName = "/GRNVBEM/ExpressionData"+str(idx)+".csv"
        outPath = 'data/' +  str(outDir) + 'outFile'+str(idx)+'.txt'

        cmdToRun = ' '.join(['docker run --rm -v', 
                             str(Path.cwd())+':/VBEM/data/ grnvbem:base /bin/sh -c \"time -v -o', 
                             "data/" + str(outDir) + 'time'+str(idx)+'.txt', 
                             './GRNVBEM', inputPath+exprName, outPath, '\"'])
        print(cmdToRun)
        os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from GRNVBEM.
    '''
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/GRNVBEM/"

    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    OutSubDF = [0]*len(colNames)

    for indx in range(len(colNames)):
        outFileName = 'outFile'+str(indx)+'.txt'
        # Quit if output file does not exist
        if not Path(outDir+outFileName).exists():
            print(outDir+outFileName+' does not exist, skipping...')
            return
        
        # Read output
        OutSubDF[indx] = pd.read_csv(outDir+outFileName, sep = '\t', header = 0)
        
    outDF = pd.concat(OutSubDF)
    FinalDF = outDF[outDF['Probability'] == outDF.groupby(['Parent','Child'])['Probability'].transform('max')]

    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for idx, row in FinalDF.sort_values(['Probability'], ascending = False).iterrows():
        outFile.write('\t'.join([row['Parent'],row['Child'],str(row['Probability'])])+'\n')
    outFile.close()
    
