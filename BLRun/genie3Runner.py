import os
import pandas as pd

from BLRun.runner import Runner


class GENIE3Runner(Runner):
    """Concrete runner for the GENIE3 GRN inference algorithm."""

    def generateInputs(self):
        '''
        Function to generate desired inputs for GENIE3.
        If the folder/files under self.input_dir exist,
        this function will not do anything.
        '''

        # Create folders in advance to prevent docker from creating folders with root-exclusive permissions
        if not self.working_dir.exists():
            self.working_dir.mkdir(parents=True, exist_ok = False)

        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True, exist_ok = False)

        # Create ExpressionData.csv file in the created input directory
        GENIE3_EXPRESSION_FILE = self.working_dir / "ExpressionData.csv"
        if not GENIE3_EXPRESSION_FILE.exists():
            # input data
            ExpressionData = pd.read_csv(self.input_dir / self.exprData,
                                         header = 0, index_col = 0)

            # Write .csv file — arboreto expects cells as rows, genes as columns
            ExpressionData.T.to_csv(GENIE3_EXPRESSION_FILE,
                                 sep = '\t', header  = True, index = True)

    def run(self):
        '''
        Function to run GENIE3 algorithm
        '''

        # Directly mount the input and output folders
        inputVolumeMount = " -v " + str(self.working_dir) + ":/input/"
        outputVolumeMount = " -v " + str(self.working_dir) + ":/output/"
        cmdToRun = ' '.join(['docker run --rm',
                            inputVolumeMount,
                            outputVolumeMount,
                            '--expose=41269',
                            'grnbeeline/arboreto:base /bin/sh -c \"time -v -o',
                            "/output/time.txt",
                            'python runArboreto.py --algo=GENIE3',
                            '--inFile=/input/ExpressionData.csv', '--outFile=/output/outFile.txt', '\"'])

        self._run_docker(cmdToRun)

    def parseOutput(self):
        '''
        Function to parse outputs from GENIE3.
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

        # Read output
        OutDF = pd.read_csv(outFile, sep = '\t', header = 0)

        self._write_ranked_edges(OutDF.rename(columns={
            'TF': 'Gene1', 'target': 'Gene2', 'importance': 'EdgeWeight'
        })[['Gene1', 'Gene2', 'EdgeWeight']])
