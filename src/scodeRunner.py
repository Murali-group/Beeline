import os
import subprocess
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

    setupParams(RunnerObj)


def setupParams(RunnerObj):
    # if the parameters aren't specified, then use default parameters
    # TODO allow passing in multiple sets of hyperparameters
    RunnerObj.params_str = build_params_str(RunnerObj.params)

    # the final file is written here:
    outDir = "outputs/" + str(RunnerObj.inputDir).split("inputs/")[1]+"/SCODE/"
    #print(outDir)
    RunnerObj.outDir = outDir
    RunnerObj.final_ranked_edges = "%s/%s-rankedEdges.tsv" % (outDir, RunnerObj.params_str)


def build_params_str(params):
    params_str = "D%s-I%s-R%s" % (params['D'], params['nIter'], params['nRep'])
    return params_str


def run(RunnerObj):
    '''
    Function to run SCODE algorithm
    '''
    
    # make output dirs if they do not exist:
    outDir = "%s/%s/" % (RunnerObj.outDir, RunnerObj.params_str)
    os.makedirs(outDir, exist_ok = True)

    inputPath = "data/"+str(RunnerObj.inputDir).split("RNMethods/")[1]+"/SCODE/"

    D = str(RunnerObj.params['D'])
    nIter = str(RunnerObj.params['nIter'])
    nRep = str(RunnerObj.params['nRep'])

    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
    colNames = PTData.columns

    for idx in range(len(colNames)):
        ExpressionData = pd.read_csv(str(RunnerObj.inputDir).split("RNMethods/")[1]+"/SCODE/"+\
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
                             nGenes, D, nCells, nIter, nRep, '\"'])
        print(cmdToRun)
        #os.system(cmdToRun)
        subprocess.check_call(cmdToRun, shell=True)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCODE.
    '''
    outDir = RunnerObj.outDir

    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
    colNames = PTData.columns
    for indx in range(len(colNames)):
        # Read output
        outFile = "%s/%s/%s/meanA.txt" % (outDir, RunnerObj.params_str, indx)
        if not Path(outFile).exists():
            # Quit if output file does not exist

            print(outFile+' does not exist, skipping...')
            return
        OutDF = pd.read_csv(outFile, sep = '\t', header = None)

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

        outFile = "%s/%s/outFile%s.csv" % (outDir, RunnerObj.params_str, indx)
        with open(outFile, 'w') as out:
            out.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')
            for row, col, val in zip(rows, cols, DFSorted):
                out.write('\t'.join([GeneList[row],GeneList[col],str(val)])+'\n')

    OutSubDF = [0]*len(colNames)
    for indx in range(len(colNames)):
        outFile = "%s/%s/outFile%s.csv" % (outDir, RunnerObj.params_str, indx)
        OutSubDF[indx] = pd.read_csv(outFile, sep = '\t', header = 0)

        OutSubDF[indx].EdgeWeight = np.abs(OutSubDF[indx].EdgeWeight)

    outDF = pd.concat(OutSubDF)
    finalDF = outDF[outDF['EdgeWeight'] == outDF.groupby(['Gene1','Gene2'])['EdgeWeight'].transform('max')]
    finalDF.sort_values(['EdgeWeight'], ascending = False, inplace = True)

    print("\twriting processed rankedEdges to %s" % (RunnerObj.final_ranked_edges))
    finalDF.to_csv(RunnerObj.final_ranked_edges, sep='\t', index=False)
