import src.scodeRunner as SCODE
import src.scnsRunner as SCNS
import src.sinceritiesRunner as SINCERITIES
import src.pidcRunner as PIDC
import src.grnvbemRunner as GRNVBEM

from pathlib import Path


InputMapper = {'SCODE':SCODE.generateInputs,
               'SINCERITIES':SINCERITIES.generateInputs,
               'SCNS':SCNS.generateInputs,
               'PIDC':PIDC.generateInputs,
               'GRNVBEM':GRNVBEM.generateInputs}

AlgorithmMapper = {'SCODE':SCODE.run,
            'SINCERITIES':SINCERITIES.run,
            'SCNS':SCNS.run,
            'PIDC':PIDC.run,
            'GRNVBEM':GRNVBEM.run}




OutputParser = {'SCODE':SCODE.parseOutput, 
            'SINCERITIES':SINCERITIES.parseOutput,
            'SCNS':SCNS.parseOutput,
            'PIDC':PIDC.parseOutput,
            'GRNVBEM':GRNVBEM.parseOutput}


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