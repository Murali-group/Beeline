import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for TENET.
    If the folder/files under RunnerObj.datadir exist,
    this function will not do anything.
    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    if not RunnerObj.inputDir.joinpath("TENET").exists():
        print("Input folder for TENET does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("TENET").mkdir(exist_ok = False)
    
    if not RunnerObj.inputDir.joinpath("SCSGL/ExpressionData.csv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath("TENET/ExpressionData.csv"),
                              sep = ',', header = True)
    
    if not RunnerObj.inputDir.joinpath("TENET/refNetwork.csv").exists():
        refNetworkData = pd.read_csv(RunnerObj.inputDir.joinpath("TENET/refNetwork.csv"),
                                     header = 0, index_col = 0)
        refNetworkData.to_csv(RunnerObj.inputDir.joinpath("TENET/refNetwork.csv"),
                              sep = ',', header = True)
    
    # TODO: Create the cell-selection raw text file for the TENET algorithm
    # if not RunnerObj.inputDir.joinpath("TENET/cell_select.txt"):
    #     pass

def run(RunnerObj):
    '''
    Function to run TENET algorithm
    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    pass

def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCSGL
    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    pass