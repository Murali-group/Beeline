import os
import pandas as pd
from pathlib import Path
import numpy as np


def generateInputs(RunnerObj):
    '''
    Function to generate desired outputs for mcpnet MCP4.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("MCP4").exists():
        print("Input folder for MCP4 does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("MCP4").mkdir(exist_ok = False)
 
    if not RunnerObj.inputDir.joinpath("MCP4/ExpressionData.csv").exists():
        import shutil
        shutil.copy(
            RunnerObj.inputDir.joinpath(RunnerObj.exprData),
            RunnerObj.inputDir.joinpath("MCP4")
        )

def run(RunnerObj):
    '''
    Function to run MCP4 algorithm
    '''
    # inputDir = str(RunnerObj.inputDir).split(str(Path.cwd()))[1]
    # inputPath = f"{inputDir}/XGBDENSE/ExpressionData.csv"
    inputPath = RunnerObj.inputDir.joinpath("MCP4/ExpressionData.csv")
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/MCP4/"
    os.makedirs(outDir, exist_ok = True)

    outPath = str(outDir) + 'outFile'
    miPath = str(outDir) + 'mi.h5'
    #TODO::
    cmdToRun = ' '.join([
        'time -v -o',
        f"{outDir}time.txt", 
        ' /bin/sh -c " mcpnet/build/bin/mi ',
        f'-i {inputPath}',
        f'-o {miPath}', 
        ';',
        'mcpnet/build/bin/mcpnet ',
        f'-i {miPath}', 
        ' --mi-input -m 1 2 3 ',
        f'-o {outPath} "', 
    ])
    print(cmdToRun)
    os.system(cmdToRun)

def parseOutput(RunnerObj):
    '''
    Function to parse outputs from MCP4.
    '''
    # Quit if output directory does not exist
    rinDir = str(RunnerObj.inputDir).split("inputs/")[1]
    outDir = f"outputs/{rinDir}/MCP4/"
    
    if not Path(outDir+'outFile_inmi.mcp4.h5').exists():
        print(outDir+'outFile_inmi.mcp4.h5'+' does not exist, skipping...')
        return
    # Read output
    h5file = outDir+'outFile_inmi.mcp4.h5'
    dfx = pd.read_hdf(h5file)
    OutDF = (
        dfx  # type:ignore
        .transpose()   # type:ignore
        .stack()
        .reset_index()
        .set_axis(
            ["TF", "target", "importance"], axis=1
        )
    )
    OutDF = OutDF[OutDF["TF"] != OutDF["target"]]
    OutDF = OutDF.sort_values(by=["importance"], ascending=False)

    final_df = OutDF.rename(columns={
        'TF': 'Gene1',
        'target': 'Gene2',
        'importance': 'EdgeWeight',
    })
    
    outPath = outDir + 'rankedEdges.csv'
    final_df.to_csv(outPath, sep='\t', index=False)
    # outFile = open(outPath,'w')
    # for idx, row in OutDF.iterrows():
    #     outFile.write('\t'.join([row['TF'],row['target'],str(row['importance'])])+'\n')
 
