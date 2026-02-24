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

        cmdToRun = ' '.join(['docker run --rm',
                            f"-v {self.working_dir}:/usr/working_dir",
                            'jump3:base /bin/sh -c \"time -v -o',
                            "/usr/working_dir/time.txt",
                            './runJump3',
                            "/usr/working_dir/ExpressionData.csv", "/usr/working_dir/outFile.txt", '\"'])

        self._run_docker(cmdToRun)

    def parseOutput(self):
        '''
        Function to parse outputs from JUMP3.
        '''
        workDir = self.working_dir
        outFile = workDir / 'outFile.txt'

        # Quit if output file does not exist
        if not outFile.exists():
            print(str(outFile) + ' does not exist, skipping...')
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

        self._write_ranked_edges(pd.DataFrame({
            'Gene1':      [GeneList[r] for r in rows],
            'Gene2':      [GeneList[c] for c in cols],
            'EdgeWeight': DFSorted,
        }))
