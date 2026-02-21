import os
import pandas as pd
import numpy as np

from BLRun.runner import Runner


class LEAPRunner(Runner):
    """Concrete runner for the LEAP GRN inference algorithm."""

    def generateInputs(self):
        '''
        Function to generate desired inputs for LEAP.
        If the folder/files under self.input_dir exist,
        this function will not do anything.
        '''

        # Create folders in advance to prevent docker from creating folders with root-exclusive permissions
        if not self.working_dir.exists():
            print("Input folder for LEAP does not exist, creating input folder...")
            self.working_dir.mkdir(parents=True, exist_ok = False)

        if not self.output_dir.exists():
            print("Output folder for LEAP does not exist, creating output folder...")
            self.output_dir.mkdir(parents=True, exist_ok = False)

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
        Function to run LEAP algorithm

        Requires the maxLag parameter
        '''

        maxLag = str(self.params['maxLag'])

        PTData = pd.read_csv(self.input_dir / self.pseudoTimeData,
                             header = 0, index_col = 0)

        colNames = PTData.columns
        for idx in range(len(colNames)):
            # Directly mount the input and output folders
            inputVolumeMount = " -v " + str(self.working_dir) + ":/input/"
            outputVolumeMount = " -v " + str(self.working_dir) + ":/output/"
            cmdToRun = ' '.join(['docker run --rm',
                                inputVolumeMount,
                                outputVolumeMount,
                                'grnbeeline/leap:base /bin/sh -c \"time -v -o',
                                "/output/time" + str(idx) + ".txt",
                                'Rscript runLeap.R',
                                "/input/ExpressionData" + str(idx) + ".csv",
                                maxLag,
                                "/output/outFile" + str(idx) + ".txt", '\"'])

            self._run_docker(cmdToRun, append=(idx > 0))

    def parseOutput(self):
        '''
        Function to parse outputs from LEAP.
        '''
        workDir = self.working_dir
        outDir = self.output_dir
        if not outDir.is_dir():
            raise FileNotFoundError(
                f"Output directory does not exist: {outDir}")

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
            OutSubDF[indx].Score = np.abs(OutSubDF[indx].Score)
        outDF = pd.concat(OutSubDF)
        FinalDF = outDF[outDF['Score'] == outDF.groupby(['Gene1','Gene2'])['Score'].transform('max')]

        self._write_ranked_edges(FinalDF.sort_values('Score', ascending=False).rename(
            columns={'Score': 'EdgeWeight'}
        )[['Gene1', 'Gene2', 'EdgeWeight']])
