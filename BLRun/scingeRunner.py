import os
import subprocess
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for SCINGE.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("SCINGE").exists():
        print("Input folder for SCINGE does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("SCINGE").mkdir(exist_ok = False)

        
        
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        # Select cells belonging to each pseudotime trajectory
        colName = colNames[idx]
        index = PTData[colName].index[PTData[colName].notnull()]
        exprName = "SCINGE/ExpressionData"+str(idx)+".csv"
        newExpressionData = ExpressionData.loc[:,index].T
        newExpressionData['PseudoTime'] = PTData.loc[index,colName]
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath(exprName),
                             sep = ',', header  = True, index = False)


def run(RunnerObj):
    '''
    Function to run SCINGE algorithm
    '''
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/SCINGE/"
    

    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SCINGE/"
    os.makedirs(outDir, exist_ok = True)

    # if the parameters aren't specified, then use default parameters
    # TODO allow passing in multiple sets of hyperparameters
    # these must be in the right order!
    params_order = [
        'lambda', 'dT', 'num_lags', 'kernel_width',
        'prob_zero_removal', 'prob_remove_samples',
        'family', 'num_replicates',
    ]
    default_params = {
        'lambda': '0.01',
        'dT': '10',
        'num_lags': '5',
        'kernel_width': '4',
        'prob_zero_removal': '0',
        'prob_remove_samples': '0.2',
        'family': 'gaussian',
        'num_replicates': '2',
    }
    params = RunnerObj.params
    for param, val in default_params.items():
        if param not in params:
            params[param] = val
    params_str = ' '.join(str(params[p]) for p in params_order) 
    
    
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):    
        outPath = str(outDir) + str(idx) + "/"
        os.makedirs(outPath, exist_ok = True)
        outFile = "data/" + outPath
        inputFile = inputPath + "ExpressionData"+str(idx)+".csv"
        
        cmdToRun = ' '.join(['docker run --rm -v', 
                             str(Path.cwd())+':/runSCINGE/data/ scinge:base /bin/sh -c \"time -v -o',
                             "data/" + str(outDir) + 'time'+str(idx)+'.txt', './runSCINGE ',
                             inputFile, outFile, params_str, '\"'])
        print(cmdToRun)
        # also print the parameters
        print("\tParameters: %s" % (', '.join("%s: %s" % (p, str(params[p])) for p in params_order)))
        subprocess.check_call(cmdToRun, shell=True)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCINGE.
    '''
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SCINGE/"
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    OutSubDF = [0]*len(colNames)

    for idx in range(len(colNames)):  
        
        # Quit if output directory does not exist
        if not Path(outDir+ str(idx)+'/SCINGE_Ranked_Edge_List.txt').exists():
            print(outDir+ str(idx)+'/SCINGE_Ranked_Edge_List.txt does not exist, skipping...')
            return

        # Read output
        OutSubDF[idx] = pd.read_csv(outDir+ str(idx)+'/SCINGE_Ranked_Edge_List.txt',
                            sep = '\t', header = 0)
    # megre the dataframe by taking the maximum value from each DF
    # Code from here: 
    # https://stackoverflow.com/questions/20383647/pandas-selecting-by-label-sometimes-return-series-sometimes-returns-dataframe
    outDF = pd.concat(OutSubDF)
    outDF.columns= ['Gene1','Gene2','EdgeWeight']
    # Group by rows code is from here:
    # https://stackoverflow.com/questions/53114609/pandas-how-to-remove-duplicate-rows-but-keep-all-rows-with-max-value
    res = outDF[outDF['EdgeWeight'] == outDF.groupby(['Gene1','Gene2'])['EdgeWeight'].transform('max')]
    # Sort values in the dataframe   
    finalDF = res.sort_values('EdgeWeight', ascending=False)   
    finalDF.to_csv(outDir+'rankedEdges.csv',sep='\t', index = False)