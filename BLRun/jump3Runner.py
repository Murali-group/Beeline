import os
import pandas as pd
import numpy as np

from BLRun.runner import Runner


class JUMP3Runner(Runner):
    """Concrete runner for the JUMP3 GRN inference algorithm."""

    def generateInputs(self):
        '''
        Function to generate desired inputs for JUMP3.
        If the folder/files under self.input_dir exist,
        this function will not do anything.
        '''

        # Create folders in advance to prevent docker from creating folders with root-exclusive permissions
        if not self.working_dir.exists():
            print("Input folder for JUMP3 does not exist, creating input folder...")
            self.working_dir.mkdir(parents=True, exist_ok = False)

        if not self.output_dir.exists():
            print("Output folder for JUMP3 does not exist, creating output folder...")
            self.output_dir.mkdir(parents=True, exist_ok = False)

        # Create ExpressionData.csv file in the created input directory
        JUMP3_EXPRESSION_FILE = self.working_dir / "ExpressionData.csv"
        if not JUMP3_EXPRESSION_FILE.exists():
            ExpressionData = pd.read_csv(self.input_dir / self.exprData,
                                         header = 0, index_col = 0)
            newExpressionData = ExpressionData.T.copy()
            PTData = pd.read_csv(self.input_dir / self.pseudoTimeData,
                                 header = 0, index_col = 0)
            # make sure the indices are strings for both dataframes
            newExpressionData.index = newExpressionData.index.map(str)
            PTData.index = PTData.index.map(str)
            # Acc. to JUMP3:
            # In input argument Time, the first time point of each time series must be 0.
            # Also has to be an integer!
            newExpressionData['Time'] = PTData['PseudoTime']-PTData['PseudoTime'].min()
            if 'Experiment' in PTData:
                newExpressionData['Experiment'] = PTData['Experiment']
            else:
                # generate it from cell number Ex_y, where x is experiment number
                #newExpressionData['Experiment'] = [int(x.split('_')[0].strip('E')) for x in PTData.index.astype(str)]
                newExpressionData['Experiment'] = 1

            newExpressionData.to_csv(JUMP3_EXPRESSION_FILE,
                                 sep = ',', header  = True, index = False)

    def run(self):
        '''
        Function to run JUMP3 algorithm
        '''

        # Directly mount the input and output folders
        inputVolumeMount = " -v " + str(self.working_dir) + ":/input/"
        outputVolumeMount = " -v " + str(self.output_dir) + ":/output/"
        cmdToRun = ' '.join(['docker run --rm',
                            inputVolumeMount,
                            outputVolumeMount,
                            'jump3:base /bin/sh -c \"time -v -o',
                            "/output/time.txt",
                            './runJump3',
                            "/input/ExpressionData.csv", "/output/outFile.txt", '\"'])

        os.system(cmdToRun)

    def parseOutput(self):
        '''
        Function to parse outputs from JUMP3.
        '''
        outDir = self.output_dir
        outFile = outDir / 'outFile.txt'
        if not outDir.is_dir():
            raise FileNotFoundError(
                f"Output directory does not exist: {outDir}")

        # Quit if output directory does not exist
        if not outFile.exists():
            print(str(outFile) +'does not exist, skipping...')
            return

        # Read output
        OutDF = pd.read_csv(outFile, sep = ',')

        # Sort values in a matrix using code from:
        # https://stackoverflow.com/questions/21922806/sort-values-of-matrix-in-python
        OutMatrix = np.abs(OutDF.values)
        idx = np.argsort(OutMatrix, axis = None)[::-1]
        rows, cols = np.unravel_index(idx, OutDF.shape)
        DFSorted = OutMatrix[rows, cols]

        # read input file for list of gene names
        ExpressionData = pd.read_csv(self.input_dir / 'ExpressionData.csv',
                                         header = 0, index_col = 0)
        GeneList = list(ExpressionData.index)

        # Write converted csv file
        outFile = open(outDir / 'rankedEdges.csv','w')
        outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

        for row, col, val in zip(rows, cols, DFSorted):
            outFile.write('\t'.join([GeneList[row],GeneList[col],str(val)])+'\n')
        outFile.close()
