import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for scSGL.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    if not RunnerObj.inputDir.joinpath("SCSGL").exists():
        print("Input folder for SCSGL does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("SCSGL").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("SCSGL/ExpressionData.csv").exists():
        # input data
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)

        # Write gene expression data in SCSGL folder 
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath("SCSGL/ExpressionData.csv"),
                             sep = ',', header  = True)

    if not RunnerObj.inputDir.joinpath("SCSGL/refNetwork.csv").exists():
        refNetworkData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges),
                                     header = 0, index_col = 0)

	  # Write reference network data in SCSGL folder 
        refNetworkData.to_csv(RunnerObj.inputDir.joinpath("SCSGL/refNetwork.csv"),
                             sep = ',', header  = True)    

    
def run(RunnerObj):
    '''
    Function to run SCSGL algorithm
    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Get path for ExpressionData.csv generated in SCSGL folder for certain type of network in inputs
    expressionDataPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/SCSGL/ExpressionData.csv"

    # Get path for refNetwor.csv generated in SCSGL folder for certain type of network in inputs
    refNetworkPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/SCSGL/refNetwork.csv"

    pos_density = str(RunnerObj.params['pos_density'])
    neg_density = str(RunnerObj.params['neg_density'])
    assoc = str(RunnerObj.params['assoc'])

    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SCSGL/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath = "data/" +  str(outDir) + 'outFile.txt'
    cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/data/ --expose=41269', 
                         'scsgl:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time.txt', 'python runScsgl.py',
                         '--expression_file='+expressionDataPath, '--ref_net_file='+refNetworkPath, '--out_file='+outPath, 
                         '--pos_density='+pos_density, '--neg_density='+neg_density, '--assoc='+assoc,
                         '\"'])

    print(cmdToRun)
    os.system(cmdToRun)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCSGL.
    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SCSGL/"

        
    # Read output file
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = '\t', header = 0)

    OutDF.sort_values(by="EdgeWeight", ascending=False, inplace=True)
    
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
    
    # Formats the outFile into a ranked edgelist comma-separated file
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for idx, row in OutDF.iterrows():
        # TODO: might need to sort
        outFile.write('\t'.join([row['Gene1'],row['Gene2'],str(row['EdgeWeight'])])+'\n')
    outFile.close()
