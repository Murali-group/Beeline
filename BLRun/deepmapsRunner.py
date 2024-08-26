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
    if not RunnerObj.inputDir.joinpath("DEEPMAPS").exists():
        print("Input folder for DEEPMAPS does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("DEEPMAPS").mkdir(exist_ok=False)

    if not RunnerObj.inputDir.joinpath("DEEPMAPS/scRNA.RDS").exists():
        # input data
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header=0, index_col=0)

        # Write .csv file
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath("DEEPMAPS/scRNA.RDS"),
                              sep='\t', header=True, index=True)

    if not RunnerObj.inputDir.joinpath("DEEPMAPS/scATAC.RDS").exists():
        # input data
        atacData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.atacData),
                               header=0, index_col=0)

        # Write .csv file
        atacData.to_csv(RunnerObj.inputDir.joinpath("DEEPMAPS/scATAC.RDS"),
                        sep='\t', header=True, index=True)

    if not RunnerObj.inputDir.joinpath("DEEPMAPS/peaks.h5").exists():
        # input data
        rna_peaks = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.rna_peaks),
                                header=0, index_col=0)

        # Write .csv file
        rna_peaks.to_csv(RunnerObj.inputDir.joinpath("DEEPMAPS/rna_peaks.h5"),
                         sep='\t', header=True, index=True)

    if not RunnerObj.inputDir.joinpath("DEEPMAPS/velo.csv").exists():
        # input data
        velo = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.velo),
                           header=0, index_col=0)

        # Write .csv file
        velo.to_csv(RunnerObj.inputDir.joinpath("DEEPMAPS/velo.csv"),
                    sep='\t', header=True, index=True)


def run(RunnerObj):
    '''
    Function to run DEEPMAPS algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + "/DEEPMAPS/"

    # make output dirs if they do not exist:
    outDir = "outputs/" + str(RunnerObj.inputDir).split("inputs/")[1] + "/DEEPMAPS/"
    os.makedirs(outDir, exist_ok=True)

    outFile = "outFile.csv"
    timeFile = "time.txt"
    velo = "velo.csv"

    if "peaks.h5":
        h5 = "rna_peaks.h5"
        cmdToRun = ' '.join(
            [
                'docker run --rm -it --init -v /Algorithms/DEEPMAPS/:/ -v /inputs/Multi_omics/DEEPMAPS:/home/user/miniconda/lib/python3.8/site-packages/lisa/data',
                '--gpus=all --ipc=host --network=host',
                str(Path.cwd()) + ':osubmbl/deepmaps-base /bin/sh -c \"time -v -o',
                "data/" + str(outDir) + timeFile,
                'Rscript --vanilla runDeepMAPS.R', '-e', inputPath + h5, '-c', inputPath + velo,
                '-o data/' + outDir, '--outFile ' + outFile])

    elif "scATAC.RDS" and "scRNA.RDS":
        exprName = "scRNA.csv"
        atacData = "scATAC.csv"
        cmdToRun = ' '.join(
            ['docker run --rm -v', str(Path.cwd()) + ':/data/grnbeeline/deepmaps:latest /bin/sh -c \"time -v -o',
             "data/" + str(outDir) + timeFile,
             'Rscript --vanilla runDeepMAPS.R', '-r', inputPath + exprName, "-a", inputPath + atacData,
             '-c', inputPath + velo, '-o data/' + outDir, '--outFile ' + outFile])

    cmdToRun += '\"'
    print(cmdToRun)
    os.system(cmdToRun)


# [docker run --rm -it --init \
#  -v /Algorithms/DEEPMAPS/:/ \
#  -v /inputs/Multi_omics/DEEPMAPS:/home/user/miniconda/lib/python3.8/site-packages/lisa/data \
#  --gpus=all \
#  --ipc=host \
#  --network=host \
#  osubmbl/deepmaps-base
#  Rscript /deepmaps/runDeepMAPS.R -  ])
def parseOutput(RunnerObj):
    '''
    Function to parse outputs from DEEPMAPS.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Quit if output directory does not exist
    outDir = "outputs/" + str(RunnerObj.inputDir).split("inputs/")[1] + "/DEEPMAPS/"

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
