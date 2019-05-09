import os
import pandas as pd
from pathlib import Path
import numpy as np


def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for SINCERITIES.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''

    setupParams(RunnerObj)

    if not RunnerObj.inputDir.joinpath("SINCERITIES").exists():
        print("Input folder for SINCERITIES does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("SINCERITIES").mkdir(exist_ok = False)
    
    
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
    colNames = PTData.columns
    for idx in range(len(colNames)):
        # Select cells belonging to each pseudotime trajectory
        colName = colNames[idx]
        index = PTData[colName].index[PTData[colName].notnull()]
        newExpressionData = ExpressionData.loc[:,index].T
        # Perform quantile binning as recommeded in the paper
        # http://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.qcut.html#pandas.qcut
        nBins = RunnerObj.params['nBins']
        tQuantiles = pd.qcut(PTData.loc[index,colName], q = nBins, duplicates ='drop')
        mid = [(a.left + a.right)/2 for a in tQuantiles]

        newExpressionData['Time'] = mid
        exprName = "SINCERITIES/nBins%s-ExpressionData%s.csv" % (nBins, idx)
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath(exprName),
                             sep = ',', header  = True, index = False)


def setupParams(RunnerObj):
    if 'nBins' not in RunnerObj.params:
        # default is 10
        RunnerObj.params['nBins'] = 10
    # TODO it only takes >= 5 
    #if RunnerObj.params['nBins'] < 5:
    #    print("ERROR")
    RunnerObj.params_str = build_params_str(RunnerObj.params)

    # the final file is written here:
    outDir = "outputs/" + str(RunnerObj.inputDir).split("inputs/")[1]+"/SINCERITIES/"
    #print(outDir)
    RunnerObj.outDir = outDir
    RunnerObj.final_ranked_edges = "%s/%s-rankedEdges.tsv" % (outDir, RunnerObj.params_str)


def build_params_str(params):
    params_str = "nBins%s" % (params['nBins'])
    return params_str

    
def run(RunnerObj):
    '''
    Function to run SINCERITIES algorithm
    '''
    inputPath = "data/" + str(RunnerObj.inputDir).split("RNMethods/")[1] + \
                        "/SINCERITIES/"
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SINCERITIES/"
    os.makedirs(outDir, exist_ok = True)
    
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        inFile = "ExpressionData"+str(idx)+".csv"
        inFile = "nBins%s-ExpressionData%s.csv" % (RunnerObj.params['nBins'], idx)
        outPath = 'data/%s/%s-outFile%s.txt' % (outDir, RunnerObj.params_str, idx)
        cmdToRun = ' '.join(['docker run --rm -v', 
                             str(Path.cwd())+':/SINCERITIES/data/ sincerities:base /bin/sh -c \"time -v -o', 
                             "data/" + str(outDir) + 'time'+str(idx)+'.txt', 'Rscript MAIN.R',
                             inputPath+inFile, outPath, '\"'])
        print(cmdToRun)
        os.system(cmdToRun)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SINCERITIES.
    '''
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SINCERITIES/"

    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
    colNames = PTData.columns
    OutSubDF = [0]*len(colNames)
    for idx in range(len(colNames)):
        # Read output
        outFile = '%s/%s-outFile%s.txt' % (outDir, RunnerObj.params_str, idx)
        if not Path(outFile).exists():
            # Quit if output file does not exist

            print(outFile+' does not exist, skipping...')
            return
        OutSubDF[idx] = pd.read_csv(outFile, sep = ',', header = 0)

    # megre the dataframe by taking the maximum value from each DF
    # From here: https://stackoverflow.com/questions/20383647/pandas-selecting-by-label-sometimes-return-series-sometimes-returns-dataframe
    outDF = pd.concat(OutSubDF)
    # Group by rows code is from here:
    # https://stackoverflow.com/questions/53114609/pandas-how-to-remove-duplicate-rows-but-keep-all-rows-with-max-value
    res = outDF[outDF['Interaction'] == outDF.groupby(['SourceGENES','TargetGENES'])['Interaction'].transform('max')]
    # Sort values in the dataframe   
    finalDF = res.sort_values('Interaction',ascending=False)
    finalDF.drop(labels = 'Edges',axis = 'columns', inplace = True)
    # SINCERITIES output is incorrectly orderd
    finalDF.columns = ['Gene2','Gene1','EdgeWeight']

    print("\twriting processed rankedEdges to %s" % (RunnerObj.final_ranked_edges))
    finalDF.to_csv(RunnerObj.final_ranked_edges,sep='\t', index = False)
    
