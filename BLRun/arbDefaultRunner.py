import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired outputs for ARBDEF.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("ARBDEF").exists():
        print("Input folder for ARBDEF does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("ARBDEF").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("ARBDEF/ExpressionData.csv").exists():
        import shutil
        shutil.copy(
            RunnerObj.inputDir.joinpath(RunnerObj.exprData),
            RunnerObj.inputDir.joinpath("ARBDEF")
        )
    
def run(RunnerObj):
    '''
    Function to run ARBDEF algorithm
    '''
    # inputDir = str(RunnerObj.inputDir).split(str(Path.cwd()))[1]
    # inputPath = f"{inputDir}/ARBDEF/ExpressionData.csv"
    inputPath = RunnerObj.inputDir.joinpath("ARBDEF/ExpressionData.csv")
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/ARBDEF/"
    os.makedirs(outDir, exist_ok = True)

    outPath = str(outDir) + 'outFile.txt'
    statsPath = str(outDir) + 'outStats.json'
    #TODO::
    cmdToRun = ' '.join([
        'time -v -o',
        f"{outDir}time.txt", 
        'python -m gbr.cli csv --method=arb:default',
        f'--out_file {outPath}', 
        f'--rstats_out_file {statsPath}',
        f'--csv_file {inputPath}',
    ])
    print(cmdToRun)
    os.system(cmdToRun)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from ARBDEF.
    '''
    # Quit if output directory does not exist
    rinDir = str(RunnerObj.inputDir).split("inputs/")[1]
    outDir = f"outputs/{rinDir}/ARBDEF/"
    
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
    # Read output
    OutDF : pd.DataFrame = pd.read_csv(
        outDir+'outFile.txt', header = 0, index_col=0
    )


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
    # outFile.close()
