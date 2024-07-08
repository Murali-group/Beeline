import os
import subprocess
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for TENET.
    If the folder/files under RunnerObj.datadir exist,
    this function will not do anything.
    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Creates input TENET directory
    if not RunnerObj.inputDir.joinpath("TENET").exists():
        print("Input folder for TENET does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("TENET").mkdir(exist_ok = False)
    
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                 header = 0, index_col = 0)
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                         header = 0, index_col = 0)
    colNames = PTData.columns
    for idx in range(len(colNames)):
        # Generates separate subdirectory for each pseudotime listing
        RunnerObj.inputDir.joinpath(f"TENET/{idx}").mkdir(exist_ok = True)
        colName = colNames[idx]
        index = PTData[colName].index[PTData[colName].notnull()]

        # File names
        exprName = f"TENET/{idx}/ExpressionData.csv"
        PTName = f"TENET/{idx}/PseudoTime.csv"
        cellSelectName = f"TENET/{idx}/CellSelection.txt"
        refNetworkName = f"TENET/{idx}/refNetwork.csv"

        # Expression Data
        ExpressionData.loc[:,index].T.to_csv(RunnerObj.inputDir.joinpath(exprName),
                                             sep = ',', header = True, index = True)
        # Pseudotime Data
        PTData.loc[index,[colName]].to_csv(str(RunnerObj.inputDir.joinpath(PTName)).split(".csv")[0] + ".txt",
                                           sep = ',', header = False, index = False)
        # Cell Selection File
        with open(RunnerObj.inputDir.joinpath(f"TENET/{idx}/CellSelection.txt"), 'w') as selectionFile:
            print(("1\n" * len(index)), end = '', file=selectionFile)
            

def run(RunnerObj):
    # # tenet wants info as gene columns and cell rows, opposite of how info is provided through pipeline
    # toTranspose = pd.read_csv(inputPath)
    # transposedDF = pd.T
    # #remove nontransposed version
    # os.remove(inputPath)
    # transposedDF.to_csv(inputPath)
    
    
    # make output dirs if they do not exist:
    
    # Get input and output main directories
    inputDir = "/data" + "/".join(str(RunnerObj.inputDir).split(str(Path.cwd()))[1].split(os.sep)) + "/TENET"
    #print(inputDir)
    outDir = "outputs/" + str(RunnerObj.inputDir).split("inputs" + os.sep)[1] + "/TENET"
    os.makedirs(outDir, exist_ok = True)

    # Get parameters
    cellHistoryLength = str(RunnerObj.params['historyLength'])
    threadsToRun = str(RunnerObj.params['threads'])

    # Get number of pseudotime listings
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData), header = 0, index_col = 0)
    colNames = PTData.columns

    for idx in range(len(colNames)):
        # Make the output directory for current Pseudotime file
        os.makedirs(f"{outDir}/{idx}/", exist_ok = True)

        # Get paths to all necessary input files
        expressionPath = f"{inputDir}/{idx}/ExpressionData.csv"
        PTPath = f"{inputDir}/{idx}/PseudoTime.txt"
        cellSelectionPath = f"{inputDir}/{idx}/CellSelection.txt"

        # MULTICORE DOES NOT WORK 
        # cmdToRun = ' '.join(
        #     [f'docker run --rm -v {Path.cwd()}:/data tenet:base', 
        #      f'/bin/sh -c "time -v -o /data/{outDir}/time{idx}.txt',
        #      f'./TENET {expressionPath} {threadsToRun}',
        #      f'{PTPath} {cellSelectionPath} {cellHistoryLength}',
        #      f'&& mv TE_result_matrix.txt /data/{outDir}/{idx}/outFile.txt"']
        # )

        cmdToRun = ' '.join(
            [f'docker run --rm -v {Path.cwd()}:/data tenet:base', 
             f'/bin/sh -c "time -v -o /data/{outDir}/time{idx}.txt',
             f'python TENETsinglecore {expressionPath} ',
             f'{PTPath} {cellSelectionPath} {cellHistoryLength}',
             f'&& mv TE_result_matrix.txt /data/{outDir}/{idx}/outFile.txt"']
        )
    
        # File path debuggers
        # print(f"Cycle: {idx}")
        # print(f"ExprPath: {expressionPath}")
        # print(f"PTPath: {PTPath}")
        # print(f"Cell Selection Path: {cellSelectionPath}")
        # print(f"Output Directory: {outDir}/{idx}/")
        # print(f"Command to run: \n{cmdToRun}\n")

        # Run command
        subprocess.check_call(cmdToRun, shell = True)

def parseOutput(RunnerObj):
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs" + os.sep)[1]+"/TENET"
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0, sep=",")
    
    colNames = PTData.columns
    OutSubDF = [0]*len(colNames)

    for idx in range(len(colNames)):
        # Read output
        outFile = f'{outDir}/{idx}/outFile.txt'
        
        if not Path(outFile).exists():
            # Quit if output file does not exist
            print(f'{outFile} does not exist, skipping...')
            return
        
        # Parse output
        parsedGRN = Path(f'{outDir}/{idx}/parsedGRN.sif')
        os.system(f'python Algorithms/TENET/TENET/makeGRNBeeline.py 1 {outFile} {parsedGRN}')
        parsedGRNFinal = Path(f'{outDir}/{idx}/parsedGRNFinal.sif')
        os.system(f'python Algorithms/TENET/TENET/trim_indirect.py {parsedGRN} 0 {parsedGRNFinal}')
        
        OutSubDF[idx] = pd.read_csv(parsedGRNFinal, sep="\t", header=None)

    # Process and finalize output
    GRN = pd.concat(OutSubDF)
    GRN.rename(columns={0: "Gene1", 1: "EdgeWeight", 2: "Gene2"}, inplace=True)
    GRN = GRN[['Gene1', 'Gene2', 'EdgeWeight']]
    GRN.sort_values('EdgeWeight', ascending=False, inplace=True)
    GRN.to_csv(f'{outDir}/rankedEdges.csv', sep="\t", index=False)