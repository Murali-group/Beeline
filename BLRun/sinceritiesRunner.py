import os
import pandas as pd

from BLRun.runner import Runner


class SINCERITIESRunner(Runner):
    """Concrete runner for the SINCERITIES GRN inference algorithm."""

    def generateInputs(self):
        '''
        Function to generate desired inputs for SINCERITIES.
        If the folder/files under self.input_dir exist,
        this function will not do anything.
        '''

        # Create folders in advance to prevent docker from creating folders with root-exclusive permissions
        if not self.working_dir.exists():
            self.working_dir.mkdir(parents=True, exist_ok = False)

        if not self.output_dir.exists():
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
            newExpressionData = ExpressionData.loc[:,index].T
            # Perform quantile binning as recommeded in the paper
            # http://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.qcut.html#pandas.qcut
            nBins = int(self.params['nBins'])
            tQuantiles = pd.qcut(PTData.loc[index,colName], q = nBins, duplicates ='drop')
            mid = [(a.left + a.right)/2 for a in tQuantiles]

            newExpressionData['Time'] = mid
            newExpressionData.to_csv(self.working_dir / exprName,
                                 sep = ',', header  = True, index = False)

    def run(self):
        '''
        Function to run SINCERITIES algorithm
        '''

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
                                'grnbeeline/sincerities:base /bin/sh -c \"time -v -o',
                                "/output/time" + str(idx) + ".txt",
                                'Rscript MAIN.R',
                                "/input/ExpressionData" + str(idx) + ".csv",
                                "/output/outFile" + str(idx) + ".txt", '\"'])

            self._run_docker(cmdToRun, append=(idx > 0))

    def parseOutput(self):
        '''
        Function to parse outputs from SINCERITIES.
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
        for idx in range(len(colNames)):
            # Read output
            outFile = 'outFile'+str(idx)+'.txt'
            if not (workDir / outFile).exists():
                # Quit if output file does not exist
                print(str(workDir / outFile) + ' does not exist, skipping...')
                return
            OutSubDF[idx] = pd.read_csv(workDir / outFile, sep = ',', header = 0)

        # megre the dataframe by taking the maximum value from each DF
        # From here: https://stackoverflow.com/questions/20383647/pandas-selecting-by-label-sometimes-return-series-sometimes-returns-dataframe
        outDF = pd.concat(OutSubDF)
        # Group by rows code is from here:
        # https://stackoverflow.com/questions/53114609/pandas-how-to-remove-duplicate-rows-but-keep-all-rows-with-max-value
        res = outDF[outDF['Interaction'] == outDF.groupby(['SourceGENES','TargetGENES'])['Interaction'].transform('max')]
        # Sort values in the dataframe
        finalDF = res.sort_values('Interaction',ascending=False)
        finalDF.drop(labels = 'Edges',axis = 'columns', inplace = True)
        # SINCERITIES output is incorrectly orderd
        finalDF.columns = ['Gene2','Gene1','EdgeWeight']
        self._write_ranked_edges(finalDF[['Gene1', 'Gene2', 'EdgeWeight']])
