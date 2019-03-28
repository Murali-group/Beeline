import os
from pathlib import Path

def generateInputs(Runner):
    '''
    Function to generate desired inputs for SCODE.
    If the folder under Runner.datadir exists, this
    function will not do anything.
    '''
    if Runner.datadir.joinpath(Runner.datasets['name']+"/SCNS").exists():
        print("Input folder for SCNS exists, skipping this step...")
    else:
        print("Input folder for SCNS does not exist, creating input folder...")

    
def run(Runner):
    '''
    Function to run SCNS algorithm
    '''
    expressionFile = 'test'
    
def parseOutput(Runner):
    '''
    Function to parse output from SCNS
    '''
    expressionFile = 'test'
    