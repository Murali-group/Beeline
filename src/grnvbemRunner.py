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
    if not RunnerObj.inputDir.joinpath("GRNVBEM").exists():
        print("Input folder for GRNVBEM does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("GRNVBEM").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("GRNVBEM/ExpressionData.csv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath('ExpressionData.csv'),
                                     header = 0, index_col = 0)
        PTData = pd.read_csv(RunnerObj.inputDir.joinpath('PseudoTime.csv'),
                             header = 0, index_col = 0)
        
        # Order columns by PseudoTime
        newExpressionData = ExpressionData[PTData.sort_values(['PseudoTime']).index]
        
        newExpressionData.insert(loc = 0, column = 'GENES', \
                                                     value = newExpressionData.index)

        # Write .csv file
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath("GRNVBEM/ExpressionData.csv"),
                             sep = ',', header  = True, index = False)
    
def run(RunnerObj):
    '''
    Function to run GRN-VBEM algorithm
    '''
    inputPath = "data/" + str(RunnerObj.inputDir).split("ModelEval/")[1] + \
                    "/GRNVBEM/ExpressionData.csv"
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/GRNVBEM/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath = "data/" +  str(outDir) + 'outFile.txt'
    cmdToRun = ' '.join(['docker run --rm -v ~/ModelEval:/VBEM/data/ grnvbem:base /bin/sh -c \"./GRNVBEM ', 
                         inputPath, outPath, '\"'])
    print(cmdToRun)
    os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from GRNVBEM.
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/GRNVBEM/"
    if not Path(outDir).exists():
        raise FileNotFoundError()
        
    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = '\t', header = 0)
    
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for idx, row in OutDF.sort_values(['Probability'], ascending = False).iterrows():
        outFile.write('\t'.join([row['Parent'],row['Child'],str(row['Probability'])])+'\n')
    outFile.close()
    