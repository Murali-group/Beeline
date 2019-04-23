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
        
    if not RunnerObj.inputDir.joinpath("SCINGE/ExpressionData.csv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        newExpressionData = ExpressionData.T.copy()
        PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
        newExpressionData['PseudoTime'] = PTData['PseudoTime']
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath("SCINGE/ExpressionData.csv"),
                             sep = ',', header  = True, index = False)


def run(RunnerObj):
    '''
    Function to run SCINGE algorithm
    '''
    inputPath = "data/" + str(RunnerObj.inputDir).split("RNMethods/")[1] + \
                    "/SCINGE/ExpressionData.csv"

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

    outPath = "data/" +  str(outDir) 
    cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/runSCINGE/data/ scinge:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time.txt', './runSCINGE ',
                         inputPath, outPath, params_str, '\"'])
    print(cmdToRun)
    # also print the parameters
    print("\tParameters: %s" % (', '.join("%s: %s" % (p, str(params[p])) for p in params_order)))
    #os.system(cmdToRun)
    subprocess.check_call(cmdToRun, shell=True)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCINGE.
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SCINGE/"
    if not Path(outDir+'SCINGE_Ranked_Edge_List.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
        
    # Read output
    OutDF = pd.read_csv(outDir+'SCINGE_Ranked_Edge_List.txt', sep = '\t', header = 0)
    
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for idx, row in OutDF.iterrows():
        outFile.write('\t'.join([row['Regulator'],row['Target'],str(row['SCINGE_Score'])])+'\n')
    outFile.close()
    
