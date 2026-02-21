import os
import pandas as pd

from BLRun.runner import Runner


class SCSGLRunner(Runner):
    """Concrete runner for the scSGL GRN inference algorithm."""

    def generateInputs(self):
        '''
        Function to generate desired inputs for scSGL.
        If the folder/files under self.input_dir exist,
        this function will not do anything.
        '''

        # Create folders in advance to prevent docker from creating folders with root-exclusive permissions
        if not self.working_dir.exists():
            print("Input folder for SCSGL does not exist, creating input folder...")
            self.working_dir.mkdir(parents=True, exist_ok = False)

        if not self.output_dir.exists():
            print("Output folder for SCSGL does not exist, creating output folder...")
            self.output_dir.mkdir(parents=True, exist_ok = False)

        # Create ExpressionData.csv file in the created input directory
        SCSGL_EXPRESSION_FILE = self.working_dir / "ExpressionData.csv"
        if not SCSGL_EXPRESSION_FILE.exists():
            # input data
            ExpressionData = pd.read_csv(self.input_dir / self.exprData,
                                         header = 0, index_col = 0)

            # Write gene expression data in SCSGL folder
            ExpressionData.to_csv(SCSGL_EXPRESSION_FILE,
                                 sep = ',', header  = True)

        SCSGL_GROUND_TRUTH_FILE = self.working_dir / "GroundTruthNetwork.csv"
        if not SCSGL_GROUND_TRUTH_FILE.exists():
            groundTruthNetworkData = pd.read_csv(self.input_dir / self.groundTruthNetwork,
                                         header = 0, index_col = 0)

            # Write reference network data in SCSGL folder
            groundTruthNetworkData.to_csv(SCSGL_GROUND_TRUTH_FILE,
                                 sep = ',', header  = True)

    def run(self):
        '''
        Function to run SCSGL algorithm
        '''

        pos_density = str(self.params['pos_density'])
        neg_density = str(self.params['neg_density'])
        assoc = str(self.params['assoc'])

        # Directly mount the input and output folders
        inputVolumeMount = " -v " + str(self.working_dir) + ":/input/"
        outputVolumeMount = " -v " + str(self.working_dir) + ":/output/"
        cmdToRun = ' '.join(['docker run --rm',
                            inputVolumeMount,
                            outputVolumeMount,
                            '--expose=41269',
                            'scsgl:base /bin/sh -c \"time -v -o',
                            "/output/time.txt", 'python run_scSGL.py',
                            '--expression_file=/input/ExpressionData.csv',
                            '--ground_truth_net_file=/input/GroundTruthNetwork.csv',
                            '--out_file=/output/outFile.txt',
                            '--pos_density='+pos_density, '--neg_density='+neg_density, '--assoc='+assoc,
                            '\"'])

        self._run_docker(cmdToRun)

    def parseOutput(self):
        '''
        Function to parse outputs from SCSGL.
        '''
        workDir = self.working_dir
        outDir = self.output_dir
        outFile = workDir / 'outFile.txt'
        if not outDir.is_dir():
            raise FileNotFoundError(
                f"Output directory does not exist: {outDir}")

        # Quit if output file does not exist
        if not outFile.exists():
            print(str(outFile) +'does not exist, skipping...')
            return

        # Read output file
        OutDF = pd.read_csv(outFile, sep = '\t', header = 0)

        OutDF.sort_values(by="EdgeWeight", ascending=False, inplace=True)

        self._write_ranked_edges(OutDF[['Gene1', 'Gene2', 'EdgeWeight']])
