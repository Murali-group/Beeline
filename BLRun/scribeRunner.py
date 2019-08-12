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
    
    
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        # Select cells belonging to each pseudotime trajectory
        colName = colNames[idx]
        index = PTData[colName].index[PTData[colName].notnull()]
        exprName = "SCRIBE/ExpressionData"+str(idx)+".csv"
        ExpressionData.loc[:,index].to_csv(RunnerObj.inputDir.joinpath(exprName),
                                 sep = ',', header  = True, index = True)
        cellName = "SCRIBE/CellData"+str(idx)+".csv"
        ptDF = PTData.loc[index,[colName]]        
        # Scribe expects a column labeled Time.
        ptDF.rename(columns = {colName:'Time'}, inplace = True)
        
        ptDF.to_csv(RunnerObj.inputDir.joinpath(cellName),
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
    
    inputPath = "data"+str(RunnerObj.inputDir).split(str(Path.cwd()))[1]+"/SCRIBE/"

    
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
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
    colNames = PTData.columns
    
    for idx in range(len(colNames)):
        # Specify file names for inputs and outputs
        exprName = "ExpressionData"+str(idx)+".csv"
        cellName = "CellData"+str(idx)+".csv"
        outFile = "outFile"+str(idx)+".csv"
        timeFile = 'time'+str(idx)+".txt"
        
        cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/data/ scribe:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + timeFile, 'Rscript runScribe.R',
                       '-e',inputPath +exprName, '-c',inputPath + cellName, 
                       '-g',inputPath + 'GeneData.csv', '-o data/'+outDir, '-d',delay, '-l', low,
                       '-m', method, '-x',fam, '--outFile '+outFile])

        if str(RunnerObj.params['log']) == 'True':
            cmdToRun += ' --log'
        if str(RunnerObj.params['ignorePT']) == 'True':
            cmdToRun += ' -i'

        cmdToRun += '\"'

        print(cmdToRun)

        os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCRIBE.
    '''
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SCRIBE/"

    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
    colNames = PTData.columns
    OutSubDF = [0]*len(colNames)
    for idx in range(len(colNames)):
        # Read output
        outFile = 'outFile'+str(idx)+'.csv'
        if not Path(outDir+outFile).exists():
            # Quit if output file does not exist

            print(outDir+outFile+' does not exist, skipping...')
            return
        OutSubDF[idx] = pd.read_csv(outDir+outFile, sep = ' ', header = None)

    # megre the dataframe by taking the maximum value from each DF
    # From here: https://stackoverflow.com/questions/20383647/pandas-selecting-by-label-sometimes-return-series-sometimes-returns-dataframe
    outDF = pd.concat(OutSubDF)
    outDF.columns= ['Gene1','Gene2','EdgeWeight']
    # Group by rows code is from here:
    # https://stackoverflow.com/questions/53114609/pandas-how-to-remove-duplicate-rows-but-keep-all-rows-with-max-value
    res = outDF[outDF['EdgeWeight'] == outDF.groupby(['Gene1','Gene2'])['EdgeWeight'].transform('max')]
    # Sort values in the dataframe   
    finalDF = res.sort_values('EdgeWeight',ascending=False)  
    
    finalDF.to_csv(outDir+'rankedEdges.csv',sep='\t', index = False)