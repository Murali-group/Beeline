import os
import subprocess
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("LEAP").exists():
        print("Input folder for LEAP does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("LEAP").mkdir(exist_ok = False)
        
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        # Select cells belonging to each pseudotime trajectory
        colName = colNames[idx]
        index = PTData[colName].index[PTData[colName].notnull()]
        exprName = "LEAP/ExpressionData"+str(idx)+".csv"
        
        subPT = PTData.loc[index,:]
        subExpr = ExpressionData[index]
        # Order columns by PseudoTime
        newExpressionData = subExpr[subPT.sort_values([colName]).index.astype(str)]
        
        newExpressionData.insert(loc = 0, column = 'GENES', \
                                                     value = newExpressionData.index)


        # Write .csv file
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath(exprName),
                             sep = ',', header  = True, index = False)

    setupParams(RunnerObj)


def setupParams(RunnerObj):
    # if the parameters aren't specified, then use default parameters
    # TODO allow passing in multiple sets of hyperparameters
    RunnerObj.params_str = build_params_str(RunnerObj.params)

    # the final file is written here:
    outDir = "outputs/" + str(RunnerObj.inputDir).split("inputs/")[1]+"/LEAP/"
    # TODO use the specified output dir
    #outDir = "%s/%s/LEAP/" % (RunnerObj.outputDir, str(RunnerObj.inputDir).split("inputs/")[1])
    #print(outDir)
    RunnerObj.outDir = outDir
    RunnerObj.final_ranked_edges = "%s/%s-rankedEdges.tsv" % (outDir, RunnerObj.params_str)


def build_params_str(params):
    params_str = "maxLag%s" % (params['maxLag'])
    return params_str


def run(RunnerObj):
    '''
    Function to run LEAP algorithm
    '''
    
    inputPath = "data/" + str(RunnerObj.inputDir).split("RNMethods/")[1]
    
    maxLag = str(RunnerObj.params['maxLag'])
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/LEAP/"
    os.makedirs(outDir, exist_ok = True)
    
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        exprName = "/LEAP/ExpressionData"+str(idx)+".csv"
        outPath = "data/%s/%s-outFile%s.txt" % (outDir, RunnerObj.params_str, idx)

        cmdToRun = ' '.join(['docker run --rm -v', 
                             str(Path.cwd())+':/data/ leap:base /bin/sh -c \"time -v -o', 
                             'data/' + str(outDir) + RunnerObj.params_str+'time'+str(idx)+'.txt', 'Rscript runLeap.R',
                             inputPath+exprName, maxLag, outPath, '\"'])
        print(cmdToRun)
        #os.system(cmdToRun)
        subprocess.check_call(cmdToRun, shell=True)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from LEAP.
    '''
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/LEAP/"

    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    OutSubDF = [0]*len(colNames)

    for indx in range(len(colNames)):
        outFile = "%s/%s-outFile%s.txt" % (outDir, RunnerObj.params_str, indx)
        # Quit if output file does not exist
        if not Path(outFile).exists():
            print(outFile+' does not exist, skipping...')
            return

        # Read output
        OutSubDF[indx] = pd.read_csv(outFile, sep = '\t', header = 0)
        # the score is the correlation.
        # In our case, a very negative correlation would indicate a high confidence edge, so take the absolute value.
        OutSubDF[indx].Score = np.abs(OutSubDF[indx].Score)
    outDF = pd.concat(OutSubDF)
    FinalDF = outDF[outDF['Score'] == outDF.groupby(['Gene1','Gene2'])['Score'].transform('max')]

    print("\twriting processed rankedEdges to %s" % (RunnerObj.final_ranked_edges))
    with open(RunnerObj.final_ranked_edges,'w') as out:
        out.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')
        for idx, row in FinalDF.sort_values(['Score'], ascending = False).iterrows():
            out.write('\t'.join([row['Gene1'],row['Gene2'],str(row['Score'])])+'\n')

