import pandas as pd
import numpy as np
def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for CNNC.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    if not RunnerObj.inputDir.joinpath("CNNC").exists():
        print("Input folder for CNNC does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("CNNC").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("CNNC/ExpressionData.csv").exists():
        # input data
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)

        # Write .csv file
        ExpressionData.T.to_csv(RunnerObj.inputDir.joinpath("CNNC/ExpressionData.csv"),
                             sep = '\t', header  = True, index = True)