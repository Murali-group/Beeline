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
    if not RunnerObj.inputDir.joinpath("PIDC").exists():
        print("Input folder for PIDC does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("PIDC").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("PIDC/ExpressionData.csv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath("PIDC/ExpressionData.csv"),
                             sep = '\t', header  = True, index = True)
    
def run(RunnerObj):
    '''
    Function to run PIDC algorithm
    '''
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/PIDC/ExpressionData.csv"
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/PIDC/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath = 'data/'+ str(outDir) + 'outFile.txt'
    cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/data pidc:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time.txt', 'julia runPIDC.jl',
                         inputPath, outPath, '\"'])
    print(cmdToRun)
    os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCODE.
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/PIDC/"
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
        
    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = '\t', header = None)
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for idx, row in OutDF.iterrows():
        outFile.write('\t'.join([row[0],row[1],str(row[2])])+'\n')
    outFile.close()
    