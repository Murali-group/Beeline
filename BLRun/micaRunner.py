import os
import pandas as pd
from pathlib import Path
import numpy as np


def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for MICA.
    If the folder/files under RunnerObj.datadir exist,
    this function will not do anything.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    if not RunnerObj.inputDir.joinpath("MICA").exists():
        print("Input folder for MICA does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("MICA").mkdir(exist_ok=False)

    if not RunnerObj.inputDir.joinpath("MICA/ExpressionData.csv").exists():
        # input data
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header=0, index_col=0)

        # Write .csv file
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath("MICA/ExpressionData.csv"),
                                sep='\t', header=True, index=True)

    if not RunnerObj.inputDir.joinpath("MICA/regulators.csv").exists():
        # input data
        reg = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.regData),
                                     header=0, index_col=0)

        # Write .csv file
        reg.to_csv(RunnerObj.inputDir.joinpath("MICA/regulators.csv"),
                                sep='\t', header=True, index=True)
    if not RunnerObj.inputDir.joinpath("MICA/atacData.csv").exists():
        # input data
        atacData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.atacData),
                                     header=0, index_col=0)

        # Write .csv file
        atacData.to_csv(RunnerObj.inputDir.joinpath("MICA/atacData.csv"),
                                sep='\t', header=True, index=True)


def run(RunnerObj):
    '''
    Function to run MICA algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + "/MICA/"

    # required inputs
    method = str(RunnerObj.params['method'])

    # make output dirs if they do not exist:
    outDir = "outputs/" + str(RunnerObj.inputDir).split("inputs/")[1] + "/MICA/"
    os.makedirs(outDir, exist_ok=True)

    exprName = "ExpressionData.csv"
    if "atacData.csv" :
        atacData = "atacData.csv"
    if "regulators.csv":
        regulators = "regulators.csv"
    outFile = "outFile.csv"
    timeFile = "time.txt"

    cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd()) + ':/data/grnbeeline/mica:latest /bin/sh -c \"time -v -o',
                         "data/" + str(outDir) + timeFile,
                         'Rscript --vanilla runMICA.R', '-r', inputPath + regulators, "-a", inputPath + atacData,
                         '-e', inputPath + exprName, '-o data/' + outDir, "-m", method, '--outFile ' + outFile])

    cmdToRun += '\"'
    print(cmdToRun)
    os.system(cmdToRun)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from MICA.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Quit if output directory does not exist
    outDir = "outputs/" + str(RunnerObj.inputDir).split("inputs/")[1] + "/MICA/"

    # Read output
    OutDF = pd.read_csv(outDir + 'outFile.csv', sep='\t', header=0)

    if not Path(outDir + 'outFile.csv').exists():
        print(outDir + 'outFile.csv' + 'does not exist, skipping...')
        return

    outFile = open(outDir + 'rankedEdges.csv', 'w')
    outFile.write('Gene1' + '\t' + 'Gene2' + '\t' + 'EdgeWeight' + '\n')

    for idx, row in OutDF.iterrows():
        outFile.write('\t'.join([row['TF'], row['target'], str(row['importance'])]) + '\n')
    outFile.close()

