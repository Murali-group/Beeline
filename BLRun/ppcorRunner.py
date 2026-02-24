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

        cmdToRun = ' '.join(['docker run --rm',
                            f"-v {self.working_dir}:/usr/working_dir",
                            'grnbeeline/ppcor:base /bin/sh -c \"time -v -o',
                            "/usr/working_dir/time.txt",
                            'Rscript runPPCOR.R',
                            "/usr/working_dir/ExpressionData.csv", "/usr/working_dir/outFile.txt", '\"'])

        # Run command
        self._run_docker(cmdToRun)

    def parseOutput(self):
        '''
        Function to parse outputs from PPCOR.
        '''
        workDir = self.working_dir
        outFile = workDir / 'outFile.txt'

        # Quit if output file does not exist
        if not outFile.exists():
            print(str(outFile) + ' does not exist, skipping...')
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
