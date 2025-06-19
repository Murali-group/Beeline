import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired outputs for XGBGPUDENSE.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("XGBGPUDENSE").exists():
        print("Input folder for XGBGPUDENSE does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("XGBGPUDENSE").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("XGBGPUDENSE/ExpressionData.csv").exists():
        import shutil
        shutil.copy(
            RunnerObj.inputDir.joinpath(RunnerObj.exprData),
            RunnerObj.inputDir.joinpath("XGBGPUDENSE")
        )
    
def run(RunnerObj):
    '''
    Function to run XGBGPUDENSE algorithm
    '''
    # inputDir = str(RunnerObj.inputDir).split(str(Path.cwd()))[1]
    # inputPath = f"{inputDir}/XGBGPUDENSE/ExpressionData.csv"
    inputPath = RunnerObj.inputDir.joinpath("XGBGPUDENSE/ExpressionData.csv")
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/XGBGPUDENSE/"
    os.makedirs(outDir, exist_ok = True)

    outPath = str(outDir) + 'outFile.txt'
    statsPath = str(outDir) + 'outStats.json'
    #TODO::
    cmdToRun = ' '.join([
        'time -v -o',
        f"{outDir}time.txt", 
        'python -m gbr.cli csv --method=xgb --device=gpu',
        f'--out_file {outPath}', 
        f'--rstats_out_file {statsPath}',
        f'--csv_file {inputPath}',
    ])
    print(cmdToRun)
    os.system(cmdToRun)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from XGBGPUDENSE.
    '''
    # Quit if output directory does not exist
    rinDir = str(RunnerObj.inputDir).split("inputs/")[1]
    outDir = f"outputs/{rinDir}/XGBGPUDENSE/"
    
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
 
