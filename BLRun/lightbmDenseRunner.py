import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired outputs for LGBDENSE.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("LGBDENSE").exists():
        print("Input folder for LGBDENSE does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("LGBDENSE").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("LGBDENSE/ExpressionData.csv").exists():
        import shutil
        shutil.copy(
            RunnerObj.inputDir.joinpath(RunnerObj.exprData),
            RunnerObj.inputDir.joinpath("LGBDENSE")
        )
    
def run(RunnerObj):
    '''
    Function to run LGBDENSE algorithm
    '''
    # inputDir = str(RunnerObj.inputDir).split(str(Path.cwd()))[1]
    # inputPath = f"{inputDir}/LGBDENSE/ExpressionData.csv"
    inputPath = RunnerObj.inputDir.joinpath("LGBDENSE/ExpressionData.csv")
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/LGBDENSE/"
    os.makedirs(outDir, exist_ok = True)

    outPath = str(outDir) + 'outFile.txt'
    statsPath = str(outDir) + 'outStats.json'
    #TODO::
    cmdToRun = ' '.join([
        'time -v -o',
        f"{outDir}time.txt", 
        'python -m gbr.cli csv --device=cpu --method=lgb',
        f'--out_file {outPath}', 
        f'--rstats_out_file {statsPath}',
        f'--csv_file {inputPath}',
    ])
    print(cmdToRun)
    os.system(cmdToRun)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from LGBDENSE.
    '''
    # Quit if output directory does not exist
    rinDir = str(RunnerObj.inputDir).split("inputs/")[1]
    outDir = f"outputs/{rinDir}/LGBDENSE/"
    
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
