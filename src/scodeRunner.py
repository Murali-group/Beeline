import os
import pandas as pd


def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for SCODE.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("SCODE").exists():
        print("Input folder for SCODE does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("SCODE").mkdir(existsok = False)
        
    if not RunnerObj.inputDir.joinpath("SCODE/ExpressionData.csv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath('ExpressionData.csv'),
                                     header = 0, index_col = 0)
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath("SCODE/ExpressionData.csv"),
                             sep = '\t', header  = False, index = False)
        
    if not RunnerObj.inputDir.joinpath("SCODE/PseudoTime.csv").exists():
        PTData = pd.read_csv(RunnerObj.inputDir.joinpath('PseudoTime.csv'),
                             header = 0, index_col = 0)
        PTData.to_csv(RunnerObj.inputDir.joinpath("SCODE/PseudoTime.csv"),
                             sep = '\t', header  = False)
        
    
def run(RunnerObj):
    '''
    Function to run SCODE algorithm
    '''
    inputPath = "data/"+str(RunnerObj.inputDir).split("ModelEval/")[1]+"/SCODE/"
    inputPath
    
    nGenes = str(RunnerObj.params['nGenes'])
    z = str(RunnerObj.params['z'])
    nCells = str(RunnerObj.params['nCells'])
    nIter = str(RunnerObj.params['nIter'])
    nRep = str(RunnerObj.params['nRep'])
    
    cmdToRun = ' '.join(['sudo docker run -v ~/ModelEval:/SCODE/data/  scode:base /bin/sh -c \"ruby run_R.rb',
                    inputPath +'/ExpressionData.csv', inputPath + '/PseudoTime.csv', nGenes, z, nCells, 
                          nIter, nRep, '\"'])
    print(cmdToRun)
    os.system(cmdToRun)
