import os

def generateInputs(Runner):
    '''
    Function to generate desired inputs for SCODE.
    If the folder under Runner.datadir exists, this
    function will not do anything.
    '''
    if Runner.inputDir.joinpath("SCODE").exists():
        print("Input folder for SCODE exists, skipping this step...")
    else:
        print("Input folder for SCODE does not exist, creating input folder...")

    
def run(Runner):
    '''
    Function to run SCODE algorithm
    '''
    expressionFile = 'test'
