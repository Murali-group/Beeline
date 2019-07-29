import os
import subprocess
import pandas as pd
from pathlib import Path
import numpy as np
from src.utils import baobab_utils
from src.scingeRunner import create_ml_lib_command


params_order = ['L', 'R', 'alphaMin']
default_params = {'L': '10', 'R': '1500', 'alphaMin': '0.3'}

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for GRISLI.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
    RunnerObj.colNames = PTData.columns

    if not RunnerObj.inputDir.joinpath("GRISLI").exists():
        print("Input folder for GRISLI does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("GRISLI").mkdir(exist_ok = False)
        
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                         header = 0, index_col = 0)
        for idx in range(len(RunnerObj.colNames)):
            RunnerObj.inputDir.joinpath("GRISLI/"+str(idx)).mkdir(exist_ok = True)
            
            # Select cells belonging to each pseudotime trajectory
            colName = RunnerObj.colNames[idx]
            index = PTData[colName].index[PTData[colName].notnull()]
            
            exprName = "GRISLI/"+str(idx)+"/ExpressionData.tsv"
            ExpressionData.loc[:,index].to_csv(RunnerObj.inputDir.joinpath(exprName),
                                     sep = '\t', header  = False, index = False)
            
            cellName = "GRISLI/"+str(idx)+"/PseudoTime.tsv"
            ptDF = PTData.loc[index,[colName]]                
            ptDF.to_csv(RunnerObj.inputDir.joinpath(cellName),
                                     sep = '\t', header  = False, index = False)

    setupParams(RunnerObj)
        

def setupParams(RunnerObj):
    # if the parameters aren't specified, then use default parameters
    # TODO allow passing in multiple sets of hyperparameters
    params = RunnerObj.params
    for param, val in default_params.items():
        if param not in params:
            params[param] = val
    params_str = build_params_str(params)
    RunnerObj.params_str = params_str
    RunnerObj.params = params

    # the final file is written here:
    outDir = "outputs/" + str(RunnerObj.inputDir).split("inputs/")[1]+"/GRISLI/"
    #print(outDir)
    RunnerObj.outDir = outDir
    RunnerObj.final_ranked_edges = "%s/%s-rankedEdges.csv" % (outDir, RunnerObj.params_str)


def build_params_str(params):
    params_str = "L%s-R%s-a%s" % (params['L'], params['R'], params['alphaMin'])
    return params_str


def run(RunnerObj):
    '''
    Function to run GRISLI algorithm
    '''

    params = RunnerObj.params
    L = str(params['L'])
    R = str(params['R'])
    alphaMin = str(params['alphaMin'])

    # if this has already been run, and forced is set to false, then skip it
    if params.get('forced') is False and \
            os.path.isfile(RunnerObj.final_ranked_edges):
        print("%s already exists. Set forced=True to overwrite" % (RunnerObj.final_ranked_edges))
        return 'already_exists'

    for idx in range(len(RunnerObj.colNames)):
        inputPath = str(RunnerObj.inputDir).split("RNMethods/")[1]+"/GRISLI/"+str(idx)+"/"
        outDir = RunnerObj.outDir+str(idx)+'/'
        # make output dirs if they do not exist:
        os.makedirs(outDir, exist_ok = True)
    
    #RunnerObj.outFile = "data/" + str(RunnerObj.outDir) + RunnerObj.params_str + '-outFile.txt'
        outFile = "%s/%s-outFile.txt" % (outDir, RunnerObj.params_str)
        #print(outFile)

        if params.get('docker') is True:
            inputPath = "data/"+inputPath
            outFile = "data/"+outFile
            cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/runGRISLI/data/ grisli:base /bin/sh -c \"time -v -o', "data/" + outDir + 'time.txt', './GRISLI ',inputPath, outFile, L, R, alphaMin,'\"'])
        else:
            grisli_path = "Algorithms/GRISLI/runGRISLI/GRISLI"
            # the time util doesn't have the -v or -o options on baobab
            # so skip them for now
            #grisli_command = "time -v -o %s %s" % (
            #    ' '.join(os.path.abspath(p) for p in [outDir+'time.txt', grisli_path, inputPath, outFile]), ' '.join([L, R, alphaMin]))
            grisli_command = "time %s %s %s %s " % (
                    grisli_path, inputPath, outFile, 
                    ' '.join([L, R, alphaMin]))
            ml_lib_command = create_ml_lib_command()
            # if this is on a cluser (e.g., baobab), write a qsub file and submit the job
            if 'qsub' in params and params['qsub'] is True:
                qsub_file = "%s/cmd.qsub" % (outDir)
                name = "grisli-%s" % (RunnerObj.params_str)
                jobs = [ml_lib_command, grisli_command]
                # TODO make the nodes, ppn and walltime parameters(?)
                baobab_utils.writeQsubFile(
                    jobs, qsub_file, name=name, nodes=1, ppn=1, walltime='10:00:00')
                cmdToRun = "qsub %s" % (qsub_file)
            else:
                cmdToRun = "%s\n%s" % (ml_lib_command, grisli_command)
        print(cmdToRun)
        #os.system(cmdToRun)
        subprocess.check_call(cmdToRun, shell=True)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from GRISLI.
    '''
    outDir = RunnerObj.outDir

    colNames = RunnerObj.colNames
    OutSubDF = [0]*len(colNames)

    for indx in range(len(colNames)):
        outFile = "%s/%s/%s-outFile.txt" % (outDir, indx, RunnerObj.params_str)
        if not Path(outFile).exists():
            print(outFile+' does not exist, skipping...')
            return
        # Read output
        OutDF = pd.read_csv(outFile, sep = ',', header = None)
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
        outFileName = "%s/%s/%s-rankedEdges.csv" % (outDir, indx, RunnerObj.params_str)
        #print("\twriting processed rankedEdges to %s" % (outFileName))
        with open(outFileName,'w') as out:
            out.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

            for row, col, val in zip(rows, cols, DFSorted):
                out.write('\t'.join([GeneList[row],GeneList[col],str((len(GeneList)*len(GeneList))-val)])+'\n')

        OutSubDF[indx] = pd.read_csv(outFileName, sep = '\t', header = 0)
    # megre the dataframe by taking the maximum value from each DF
    # From here: https://stackoverflow.com/questions/20383647/pandas-selecting-by-label-sometimes-return-series-sometimes-returns-dataframe
    outDF = pd.concat(OutSubDF)

    res = outDF.groupby(['Gene1','Gene2'],as_index=False).mean()
    # Sort values in the dataframe   
    finalDF = res.sort_values('EdgeWeight',ascending=False)  
    
    print("\twriting processed rankedEdges to %s" % (RunnerObj.final_ranked_edges))
    finalDF.to_csv(RunnerObj.final_ranked_edges,sep='\t', index = False)
