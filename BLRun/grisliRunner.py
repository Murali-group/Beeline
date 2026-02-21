import os
import pandas as pd
import numpy as np

from BLRun.runner import Runner


class GRISLIRunner(Runner):
    """Concrete runner for the GRISLI GRN inference algorithm."""

    def generateInputs(self):
        '''
        Function to generate desired inputs for GRISLI.
        If the folder/files under self.input_dir exist,
        this function will not do anything.
        '''

        # Create folders in advance to prevent docker from creating folders with root-exclusive permissions
        if not self.working_dir.exists():
            print("Input folder for GRISLI does not exist, creating input folder...")
            self.working_dir.mkdir(parents=True, exist_ok = False)

        if not self.output_dir.exists():
            print("Output folder for GRISLI does not exist, creating output folder...")
            self.output_dir.mkdir(parents=True, exist_ok = False)

        ExpressionData = pd.read_csv(self.input_dir / self.exprData,
                                         header = 0, index_col = 0)
        PTData = pd.read_csv(self.input_dir / self.pseudoTimeData,
                             header = 0, index_col = 0)

        colNames = PTData.columns
        for idx in range(len(colNames)):
            (self.working_dir / str(idx)).mkdir(exist_ok = True)

            # Select cells belonging to each pseudotime trajectory
            colName = colNames[idx]
            index = PTData[colName].index[PTData[colName].notnull()]

            exprName = str(idx)+"/ExpressionData.tsv"
            ExpressionData.loc[:,index].to_csv(self.working_dir / exprName,
                                     sep = '\t', header  = False, index = False)

            cellName = str(idx)+"/PseudoTime.tsv"
            ptDF = PTData.loc[index,[colName]]
            ptDF.to_csv(self.working_dir / cellName,
                                     sep = '\t', header  = False, index = False)

    def run(self):
        '''
        Function to run GRISLI algorithm
        '''

        L = str(self.params['L'])
        R = str(self.params['R'])
        alphaMin = str(self.params['alphaMin'])

        PTData = pd.read_csv(self.input_dir / self.pseudoTimeData,
                             header = 0, index_col = 0)

        colNames = PTData.columns
        for idx in range(len(colNames)):
            os.makedirs(str(self.working_dir / str(idx)), exist_ok = True)

            # Directly mount the input and output folders
            inputVolumeMount = " -v " + str(self.working_dir) + ":/input/"
            outputVolumeMount = " -v " + str(self.working_dir) + ":/output/"
            cmdToRun = ' '.join(['docker run --rm',
                                inputVolumeMount,
                                outputVolumeMount,
                                'grnbeeline/grisli:base /bin/sh -c \"time -v -o',
                                "/output/time" + str(idx) + ".txt",
                                './GRISLI',
                                "/input/" + str(idx) + "/",
                                "/output/" + str(idx) + "/outFile.txt",
                                L, R, alphaMin, '\"'])

            os.system(cmdToRun)

    def parseOutput(self):
        '''
        Function to parse outputs from GRISLI.
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
            # Read output
            outFile = str(indx)+'/outFile.txt'
            if not (workDir / outFile).exists():
                # Quit if output file does not exist
                print(str(workDir / outFile) + ' does not exist, skipping...')
                return
            OutDF = pd.read_csv(workDir / outFile, sep = ',', header = None)
            # Sort values in a matrix using code from:
            # https://stackoverflow.com/questions/21922806/sort-values-of-matrix-in-python
            OutMatrix = OutDF.values
            idx = np.argsort(OutMatrix, axis = None)
            rows, cols = np.unravel_index(idx, OutDF.shape)
            DFSorted = OutMatrix[rows, cols]

            # read input file for list of gene names
            ExpressionData = pd.read_csv(self.input_dir / self.exprData,
                                             header = 0, index_col = 0)
            GeneList = list(ExpressionData.index)
            outFileName = workDir / str(indx) / 'rankedEdges.csv'
            outFile = open(outFileName,'w')
            outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

            for row, col, val in zip(rows, cols, DFSorted):
                outFile.write('\t'.join([GeneList[row],GeneList[col],str((len(GeneList)*len(GeneList))-val)])+'\n')
            outFile.close()

            OutSubDF[indx] = pd.read_csv(outFileName, sep = '\t', header = 0)

            # megre the dataframe by taking the maximum value from each DF
            # From here: https://stackoverflow.com/questions/20383647/pandas-selecting-by-label-sometimes-return-series-sometimes-returns-dataframe
        outDF = pd.concat(OutSubDF)

        res = outDF.groupby(['Gene1','Gene2'],as_index=False).max()
        # Sort values in the dataframe
        finalDF = res.sort_values('EdgeWeight',ascending=False)

        self._write_ranked_edges(finalDF)
