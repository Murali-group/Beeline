import os
import pandas as pd

from BLRun.runner import Runner


class GRNVBEMRunner(Runner):
    """Concrete runner for the GRN-VBEM GRN inference algorithm."""

    def generateInputs(self):
        '''
        Function to generate desired inputs for GRNVBEM.
        It will create the input folder at self.working_dir if it
        does not exist already. The input folder will contain an ExpressionData.csv with
        cells ordered according to the pseudotime along the columns, and genes along
        the rows. If the files already exist, this function will overwrite it.
        '''

        ExpressionData = pd.read_csv(self.input_dir / self.exprData,
                                         header = 0, index_col = 0)
        PTData = pd.read_csv(self.input_dir / self.pseudoTimeData,
                             header = 0, index_col = 0)

        colNames = PTData.columns
        for idx in range(len(colNames)):
            # Select cells belonging to each pseudotime trajectory
            colName = colNames[idx]
            index = PTData[colName].index[PTData[colName].notnull()]
            exprName = "ExpressionData"+str(idx)+".csv"

            subPT = PTData.loc[index,:]
            subExpr = ExpressionData[index]
            # Order columns by PseudoTime
            newExpressionData = subExpr[subPT.sort_values([colName]).index.astype(str)]

            newExpressionData.insert(loc = 0, column = 'GENES', \
                                                         value = newExpressionData.index)

            # Write .csv file
            newExpressionData.to_csv(self.working_dir / exprName,
                                 sep = ',', header  = True, index = False)

    def run(self):
        '''
        Function to run GRN-VBEM algorithm
        '''

        PTData = pd.read_csv(self.input_dir / self.pseudoTimeData,
                             header = 0, index_col = 0)

        colNames = PTData.columns
        for idx in range(len(colNames)):
            cmdToRun = ' '.join(['docker run --rm',
                                f"-v {self.working_dir}:/usr/working_dir",
                                'grnbeeline/grnvbem:base /bin/sh -c \"time -v -o',
                                "/usr/working_dir/time" + str(idx) + ".txt",
                                './GRNVBEM',
                                "/usr/working_dir/ExpressionData" + str(idx) + ".csv",
                                "/usr/working_dir/outFile" + str(idx) + ".txt", '\"'])

            self._run_docker(cmdToRun, append=(idx > 0))

    def parseOutput(self):
        '''
        Function to parse outputs from GRNVBEM.
        '''
        workDir = self.working_dir

        PTData = pd.read_csv(self.input_dir / self.pseudoTimeData,
                             header = 0, index_col = 0)

        colNames = PTData.columns
        OutSubDF = [0]*len(colNames)

        for indx in range(len(colNames)):
            outFileName = 'outFile'+str(indx)+'.txt'
            # Quit if output file does not exist
            if not (workDir / outFileName).exists():
                print(str(workDir / outFileName) + ' does not exist, skipping...')
                return

            # Read output
            OutSubDF[indx] = pd.read_csv(workDir / outFileName, sep = '\t', header = 0)

        outDF = pd.concat(OutSubDF)
        FinalDF = outDF[outDF['Probability'] == outDF.groupby(['Parent','Child'])['Probability'].transform('max')]

        self._write_ranked_edges(
            FinalDF.sort_values('Probability', ascending=False).rename(
                columns={'Parent': 'Gene1', 'Child': 'Gene2', 'Probability': 'EdgeWeight'}
            )[['Gene1', 'Gene2', 'EdgeWeight']]
        )
