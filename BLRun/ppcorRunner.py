import os
import pandas as pd

from BLRun.runner import Runner


class PPCORRunner(Runner):
    """Concrete runner for the PPCOR GRN inference algorithm."""

    def generateInputs(self):
        '''
        Function to generate desired inputs for PPCOR.
        If the folder/files under self.input_dir exist,
        this function will not do anything.
        '''

        # Create folders in advance to prevent docker from creating folders with root-exclusive permissions
        if not self.working_dir.exists():
            print("Input folder for PPCOR does not exist, creating input folder...")
            self.working_dir.mkdir(parents=True, exist_ok = False)

        if not self.output_dir.exists():
            print("Output folder for PPCOR does not exist, creating output folder...")
            self.output_dir.mkdir(parents=True, exist_ok = False)

        # Create ExpressionData.csv file in the created input directory
        PPCOR_EXPRESSION_FILE = self.working_dir / "ExpressionData.csv"
        if not PPCOR_EXPRESSION_FILE.exists():
            ExpressionData = pd.read_csv(self.input_dir / self.exprData,
                                         header = 0, index_col = 0)

            newExpressionData = ExpressionData.copy()

            # Write .csv file
            newExpressionData.to_csv(PPCOR_EXPRESSION_FILE,
                                 sep = ',', header  = True, index = True)

    def run(self):
        '''
        Function to run PPCOR algorithm
        '''

        # Directly mount the input and output folders
        inputVolumeMount = " -v " + str(self.working_dir) + ":/input/"
        outputVolumeMount = " -v " + str(self.working_dir) + ":/output/"
        cmdToRun = ' '.join(['docker run --rm',
                            inputVolumeMount,
                            outputVolumeMount,
                            'grnbeeline/ppcor:base /bin/sh -c \"time -v -o',
                            "/output/time.txt",
                            'Rscript runPPCOR.R',
                            "/input/ExpressionData.csv", "/output/outFile.txt", '\"'])

        # Run command
        os.system(cmdToRun)

    def parseOutput(self):
        '''
        Function to parse outputs from PPCOR.
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
        # edges with significant p-value
        part1 = OutDF.loc[OutDF['pValue'] <= float(self.params['pVal'])]
        part1 = part1.assign(absCorVal = part1['corVal'].abs())
        # edges without significant p-value
        part2 = OutDF.loc[OutDF['pValue'] > float(self.params['pVal'])]

        part1_sorted = part1.sort_values('absCorVal', ascending=False)
        part2_out = part2[['Gene1', 'Gene2']].copy()
        part2_out['EdgeWeight'] = 0.0

        self._write_ranked_edges(pd.concat([
            part1_sorted[['Gene1', 'Gene2']].assign(EdgeWeight=part1_sorted['corVal']),
            part2_out,
        ], ignore_index=True)[['Gene1', 'Gene2', 'EdgeWeight']])
