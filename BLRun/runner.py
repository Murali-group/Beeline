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
import os
import subprocess
import shutil
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
               'SCSGL':SCSGL.generateInputs}




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
            'SCSGL':SCSGL.run}



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
            'SCSGL':SCSGL.parseOutput}


class Runner(object):
    '''
    A runnable analysis to be incorporated into the pipeline
    '''
    def __init__(self, params):
        self.name = params['name']
        self.inputDir = params['inputDir']
        self.params = params['params']
        self.exprData = params['exprData']
        self.cellData = params['cellData']
        self.trueEdges = params['trueEdges']  # used for evaluation
        self.use_embeddings = params.get('use_embeddings', False)
        self.embeddings_file = None
        
    def generate_embeddings(self):
        embed_script_path = Path(__file__).resolve().parent / "generate_embeds.py"
        
        if not embed_script_path.exists():
            raise FileNotFoundError(f"Embeddings script not found at {embed_script_path}")
        
        expr_filename = Path(self.exprData).stem
        new_input_dir = Path(self.inputDir) / f"processed_{expr_filename}"
        new_input_dir.mkdir(parents=True, exist_ok=True)
        
        self.embeddings_file = new_input_dir / "EmbeddingsData.csv"
        
        print(f"Using input file: {self.exprData}")
        print(f"Running embedding generation script at: {embed_script_path}")
        print(f"Embeddings will be saved to: {self.embeddings_file}")
        
        command = [
            "python",
            str(embed_script_path),
            "--input", str(Path(self.inputDir) / self.exprData),
        ]
        
        try:
            subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print("Embeddings generated successfully.")
            
            generated_file = Path(self.inputDir) / "EmbeddingsData.csv"
            if generated_file.exists():
                shutil.move(str(generated_file), str(self.embeddings_file))
                print(f"Embeddings moved to {self.embeddings_file}")
            
            ref_network_file = Path(self.inputDir) / "refNetwork.csv"
            if ref_network_file.exists():
                shutil.copy(str(ref_network_file), str(new_input_dir))
                print(f"refNetwork.csv copied to {new_input_dir}")
            else:
                print("refNetwork.csv not found in the input directory.")
            
            self.inputDir = new_input_dir
            self.exprData = "EmbeddingsData.csv"
            
        except subprocess.CalledProcessError as e:
            print(f"Error generating embeddings: {e}")
            print(f"Script output: {e.output.decode()}")
            print(f"Script error: {e.stderr.decode()}")
            raise

    def generateInputs(self):
        if self.use_embeddings:
            self.generate_embeddings()
        InputMapper[self.name](self)
        
    def run(self):
        input_path = Path(self.inputDir)
        output_path = Path(str(input_path).replace("inputs", "outputs"))
        output_path.mkdir(parents=True, exist_ok=True)
        
        AlgorithmMapper[self.name](self)

    def parseOutput(self):
        OutputParser[self.name](self)