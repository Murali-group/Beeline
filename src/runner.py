import src.scodeRunner as SCODE
import src.scnsRunner as SCNS
from pathlib import Path


InputMapper = {'SCODE':SCODE.generateInputs,
            'SCNS':SCNS.generateInputs}

AlgorithmMapper = {'SCODE':SCODE.run,
            'SCNS':SCNS.run}


OutputParser = {'SCODE':SCODE.parseOutput, 
            'SCNS':SCNS.parseOutput}

class Runner(object):
    '''
    A runnable analysis to be incorporated into the pipeline
    '''
    def __init__(self,
                params):
        self.name = params['name']
        self.inputDir = params['inputDir']
        self.params = params['params']
        
    def generateInputs(self):
        InputMapper[self.name](self)
        
        
    def run(self):
        AlgorithmMapper[self.name](self)


    def parseOutput(self):
        OutputParser[self.name](self)