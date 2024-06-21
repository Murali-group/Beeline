import os
import pandas as pd
from pathlib import Path
import numpy as np
import csv
def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for TENET.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("TENET").exists():
        print("Input folder for TENET does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("TENET").mkdir(exist_ok = False)


    inputPath = "data"+str(RunnerObj.inputDir).split(str(Path.cwd()))[1]+"/TENET/"
    os.makedirs(inputPath,exist_ok = True ) 

    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData)
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    #TODO write the pseudotime value in PTData one line at a time for each cell to an output trajectory.txt file 
    filename = 'trajectory.txt'

    path = inputPath + 'trajectory.txt'
    
    #export DataFrame to text file
    with open(path+'trajectory.txt', 'a') as f:
       df_string = PTData.to_string(header=False, index=False)
       f.write(df_string)
    
   #TODO write to a expression_data.txt file
    filename2 = "expression_data.txt" 
    ExpressionData = ExpressionData.transpose()
    path = inputPath + 'expression_data.txt'
                                 
    #export DataFrame to text file
    with open(path, 'a') as f:
        df_string = ExpressionData.to_string(header=False, index=False)
        f.write(df_string)

    #TODO write a cell_select.txt file with a 1 in each line for each cell
    #either the rows/columns for the ExpressionData
   
    with open(inputPath + 'cell.select.txt', 'w') as thirdfile:
        for line in range(ExpressionData.shape[0]):
            thirdfile.write('1')
    

    #Convert BEELINE ExpressionData.csv and Pseudotime.csv to the file format expected by TENET
    #which looks like https://github.com/neocaleb/TENET/blob/master/expression_data.csv and
    #https://github.com/neocaleb/TENET/blob/master/pseudotimeTuck.txt
    
def run(RunnerObj):
    #Function to run TENET algorithm
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/TENET/"
    os.makedirs(outDir, exist_ok = True)
    
    inputPath = "data"+str(RunnerObj.inputDir).split(str(Path.cwd()))[1]+"/TENET/"

    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    ExpressionData = pd.read_csv(str(RunnerObj.inputDir).split(str(Path.cwd())+'/')[1]+"/TENET/"+\
                                     "ExpressionData.csv",
                                     header = None, index_col = None, sep ='\t')
    nCells = str(ExpressionData.shape[1])
    nGenes = str(ExpressionData.shape[0])

    os.makedirs(outDir, exist_ok = True)
    
    cmdToRun = ' '.join(['docker run --rm --workdir /TENET grnbeeline/tenet:base', inputPath + 'expressionData.txt', '10', inputPath + 'trajectory.txt',inputPath +  'cell_select.txt',$
    print(cmdToRun)
    os.system(cmdToRun)

def parseOutput(RunnerObj):
    '''
    #Function to parse outputs from TENET.
    '''
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/TENET/"

    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
    colNames = PTData.columns
    # Read output

    #TODO read output from TENET and convert to a pandas dataframe DFSorted

    # read input file for list of gene names
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
    GeneList = list(ExpressionData.index)

    #Got to here
                                 
