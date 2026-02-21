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

        # Create folders in advance to prevent docker from creating folders with root-exclusive permissions
        if not self.working_dir.exists():
            print("Input folder for GRNVBEM does not exist, creating input folder...")
            self.working_dir.mkdir(parents=True, exist_ok = False)

        if not self.output_dir.exists():
            print("Output folder for GRNVBEM does not exist, creating output folder...")
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
        Function to run GRN-VBEM algorithm
        '''

        PTData = pd.read_csv(self.input_dir / self.pseudoTimeData,
                             header = 0, index_col = 0)

        colNames = PTData.columns
        for idx in range(len(colNames)):
            # Directly mount the input and output folders
            inputVolumeMount = " -v " + str(self.working_dir) + ":/input/"
            outputVolumeMount = " -v " + str(self.output_dir) + ":/output/"
            cmdToRun = ' '.join(['docker run --rm',
                                inputVolumeMount,
                                outputVolumeMount,
                                'grnbeeline/grnvbem:base /bin/sh -c \"time -v -o',
                                "/output/time" + str(idx) + ".txt",
                                './GRNVBEM',
                                "/input/ExpressionData" + str(idx) + ".csv",
                                "/output/outFile" + str(idx) + ".txt", '\"'])

            os.system(cmdToRun)

    def parseOutput(self):
        '''
        Function to parse outputs from GRNVBEM.
        '''
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
            if not (outDir / outFileName).exists():
                print(str(outDir / outFileName) + ' does not exist, skipping...')
                return

            # Read output
            OutSubDF[indx] = pd.read_csv(outDir / outFileName, sep = '\t', header = 0)

        outDF = pd.concat(OutSubDF)
        FinalDF = outDF[outDF['Probability'] == outDF.groupby(['Parent','Child'])['Probability'].transform('max')]

        # Write converted csv file
        outFile = open(outDir / 'rankedEdges.csv','w')
        outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

        for _, row in FinalDF.sort_values(['Probability'], ascending = False).iterrows():
            outFile.write('\t'.join([row['Parent'],row['Child'],str(row['Probability'])])+'\n')
        outFile.close()
