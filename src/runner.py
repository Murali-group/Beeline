import src.scodeRunner as SCODE
import src.scnsRunner as SCNS
import src.sinceritiesRunner as SINCERITIES
import src.pidcRunner as PIDC
import src.grnvbemRunner as GRNVBEM
import src.genie3Runner as GENIE3
import src.grnboost2Runner as GRNBOOST2
import src.leapRunner as LEAP
import src.jump3Runner as JUMP3
import src.ppcorRunner as PPCOR
import src.grisliRunner as GRISLI
import src.scingeRunner as SCINGE
import src.scribeRunner as SCRIBE

from pathlib import Path


LibMapper = {
    'SCODE':SCODE,
    'SINCERITIES':SINCERITIES,
    'SCNS':SCNS,
    'PIDC':PIDC,
    'GRNVBEM':GRNVBEM,
    'GENIE3':GENIE3,
    'GRNBOOST2':GRNBOOST2,
    'LEAP':LEAP,
    'JUMP3':JUMP3,
    'PPCOR':PPCOR,
    'GRISLI':GRISLI,
    'SCINGE':SCINGE,
    'SCRIBE':SCRIBE}



class Runner(object):
    '''
    A runnable analysis to be incorporated into the pipeline
    '''
    def __init__(self,
                params):
        self.name = params['name']
        self.inputDir = params['inputDir']
        self.outputDir = params['outputDir']
        self.params = params['params']
        self.exprData = params['exprData']
        self.cellData = params['cellData']

    def generateInputs(self):
        LibMapper[self.name].generateInputs(self)

    def run(self):
        return LibMapper[self.name].run(self)

    def parseOutput(self):
        LibMapper[self.name].parseOutput(self)

    def setupParams(self):
        LibMapper[self.name].setupParams(self)
