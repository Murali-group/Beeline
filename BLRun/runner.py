import BLRun.scodeRunner as SCODE
import BLRun.scnsRunner as SCNS
import BLRun.sinceritiesRunner as SINCERITIES
import BLRun.pidcRunner as PIDC
import BLRun.grnvbemRunner as GRNVBEM
import BLRun.genie3Runner as GENIE3
import BLRun.grnboost2Runner as GRNBOOST2
import BLRun.leapRunner as LEAP
import BLRun.jump3Runner as JUMP3
import BLRun.ppcorRunner as PPCOR
import BLRun.grisliRunner as GRISLI
import BLRun.singeRunner as SINGE
import BLRun.scribeRunner as SCRIBE
import BLRun.scsglRunner as SCSGL
import BLRun.xgbDenseRunner as XGBDR
import BLRun.xgbGPUDenseRunner as XGBGPDR
import BLRun.lightbmDenseRunner as LBMDR
import BLRun.arbDefaultRunner as ARBDR
import BLRun.miRunner as MIRR
import BLRun.clrRunner as CLRR
import BLRun.dpiRunner as DPIR
import BLRun.mcp2Runner as MCP2R
import BLRun.mcp3Runner as MCP3R
import BLRun.mcp4Runner as MCP4R

from pathlib import Path

InputMapper = {'SCODE':SCODE.generateInputs,
               'SINCERITIES':SINCERITIES.generateInputs,
               'SCNS':SCNS.generateInputs,
               'PIDC':PIDC.generateInputs,
               'GRNVBEM':GRNVBEM.generateInputs,
               'GENIE3':GENIE3.generateInputs,
               'GRNBOOST2':GRNBOOST2.generateInputs,
               'LEAP':LEAP.generateInputs,
               'JUMP3':JUMP3.generateInputs,
               'PPCOR':PPCOR.generateInputs,
               'GRISLI':GRISLI.generateInputs,
               'SINGE':SINGE.generateInputs,
               'SCRIBE':SCRIBE.generateInputs,
               'SCSGL':SCSGL.generateInputs,
               'XGBDENSE': XGBDR.generateInputs,
               'XGBGPUDENSE': XGBGPDR.generateInputs,
               'LGBDENSE': LBMDR.generateInputs,
               'ARBDEF': ARBDR.generateInputs,
               'MI': MIRR.generateInputs,
               'CLR': CLRR.generateInputs,
               'DPI': DPIR.generateInputs,
               'MCP2': MCP2R.generateInputs,
               'MCP3': MCP3R.generateInputs,
               'MCP4': MCP4R.generateInputs,
               }




AlgorithmMapper = {'SCODE':SCODE.run,
            'SINCERITIES':SINCERITIES.run,
            'SCNS':SCNS.run,
            'PIDC':PIDC.run,
            'GRNVBEM':GRNVBEM.run,
            'GENIE3':GENIE3.run,
            'GRNBOOST2':GRNBOOST2.run,
            'LEAP':LEAP.run,
            'JUMP3':JUMP3.run,
            'PPCOR':PPCOR.run,
            'GRISLI':GRISLI.run,
            'SINGE':SINGE.run,
            'SCRIBE':SCRIBE.run,
            'SCSGL':SCSGL.run,
            'XGBDENSE':XGBDR.run,
            'XGBGPUDENSE':XGBGPDR.run,
            'LGBDENSE': LBMDR.run,
            'ARBDEF':ARBDR.run,
            'MI': MIRR.run,
            'CLR': CLRR.run,
            'DPI': DPIR.run,
            'MCP2': MCP2R.run,
            'MCP3': MCP3R.run,
            'MCP4': MCP4R.run,
            }



OutputParser = {'SCODE':SCODE.parseOutput, 
            'SINCERITIES':SINCERITIES.parseOutput,
            'SCNS':SCNS.parseOutput,
            'PIDC':PIDC.parseOutput,
            'GRNVBEM':GRNVBEM.parseOutput,
            'GENIE3':GENIE3.parseOutput,
            'GRNBOOST2':GRNBOOST2.parseOutput,
            'LEAP': LEAP.parseOutput,
            'JUMP3': JUMP3.parseOutput,
            'PPCOR':PPCOR.parseOutput,
            'GRISLI':GRISLI.parseOutput,
            'SINGE':SINGE.parseOutput,
            'SCRIBE':SCRIBE.parseOutput,
            'SCSGL':SCSGL.parseOutput,
            'XGBDENSE':XGBDR.parseOutput,
            'XGBGPUDENSE':XGBGPDR.parseOutput,
            'LGBDENSE': LBMDR.parseOutput,
            'ARBDEF':ARBDR.parseOutput,
            'MI': MIRR.parseOutput,
            'CLR': CLRR.parseOutput,
            'DPI': DPIR.parseOutput,
            'MCP2': MCP2R.parseOutput,
            'MCP3': MCP3R.parseOutput,
            'MCP4': MCP4R.parseOutput,
            }


class Runner(object):
    '''
    A runnable analysis to be incorporated into the pipeline
    '''
    def __init__(self,
                params):
        self.name = params['name']
        self.inputDir = params['inputDir']
        self.params = params['params']
        self.exprData = params['exprData']
        self.cellData = params['cellData']
        self.trueEdges = params['trueEdges'] #used for evaluation
        
    def generateInputs(self):
        InputMapper[self.name](self)
        
        
    def run(self):
        AlgorithmMapper[self.name](self)

    def parseOutput(self):
        OutputParser[self.name](self)
