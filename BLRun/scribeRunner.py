import os
import pandas as pd

from BLRun.runner import Runner


class SCRIBERunner(Runner):
    """Concrete runner for the SCRIBE GRN inference algorithm."""

    def generateInputs(self):
        '''
        Function to generate desired inputs for SCRIBE.
        If the folder/files under self.input_dir exist,
        this function will not do anything.
        '''

        # Create folders in advance to prevent docker from creating folders with root-exclusive permissions
        if not self.working_dir.exists():
            print("Input folder for SCRIBE does not exist, creating input folder...")
            self.working_dir.mkdir(parents=True, exist_ok = False)

        if not self.output_dir.exists():
            print("Output folder for SCRIBE does not exist, creating output folder...")
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
            ExpressionData.loc[:,index].to_csv(self.working_dir / exprName,
                                     sep = ',', header  = True, index = True)
            cellName = "pseudoTimeData"+str(idx)+".csv"
            ptDF = PTData.loc[index,[colName]]
            # Scribe expects a column labeled Time.
            ptDF.rename(columns = {colName:'Time'}, inplace = True)

            ptDF.to_csv(self.working_dir / cellName,
                                     sep = ',', header  = True, index = True)

        SCRIBE_GENE_FILE = self.working_dir / "GeneData.csv"
        if not SCRIBE_GENE_FILE.exists():
            # required column!!
            geneDict = {}
            geneDict['gene_short_name'] = [gene.replace('x_', '') for gene in ExpressionData.index]

            geneDF = pd.DataFrame(geneDict, index = ExpressionData.index)
            geneDF.to_csv(SCRIBE_GENE_FILE,
                          sep = ',', header = True)

    def run(self):
        '''
        Function to run SCRIBE algorithm.
        To see all the inputs runScribe.R script takes, run:
        docker run scribe:base /bin/sh -c "Rscript runScribe.R -h"
        '''

        # required inputs
        delay = str(self.params['delay'])
        method = str(self.params['method'])
        low = str(self.params['lowerDetectionLimit'])
        fam = str(self.params['expressionFamily'])

        # Build the command to run Scribe
        PTData = pd.read_csv(self.input_dir / self.pseudoTimeData,
                             header = 0, index_col = 0)
        colNames = PTData.columns

        for idx in range(len(colNames)):
            # Specify file names for inputs and outputs
            exprName = "ExpressionData"+str(idx)+".csv"
            cellName = "pseudoTimeData"+str(idx)+".csv"
            outFile = "outFile"+str(idx)+".csv"
            timeFile = 'time'+str(idx)+".txt"

            # Directly mount the input and output folders
            inputVolumeMount = " -v " + str(self.working_dir) + ":/input/"
            outputVolumeMount = " -v " + str(self.working_dir) + ":/output/"
            cmdToRun = ' '.join(['docker run --rm',
                           inputVolumeMount,
                           outputVolumeMount,
                           'grnbeeline/scribe:base /bin/sh -c \"time -v -o',
                           "/output/" + timeFile, 'Rscript runScribe.R',
                           '-e', "/input/" + exprName, '-c', "/input/" + cellName,
                           '-g', "/input/GeneData.csv", '-o /output/', '-d', delay, '-l', low,
                           '-m', method, '-x', fam, '--outFile ' + outFile])

            if str(self.params['log']) == 'True':
                cmdToRun += ' --log'
            if str(self.params['ignorePT']) == 'True':
                cmdToRun += ' -i'

            cmdToRun += '\"'

            self._run_docker(cmdToRun, append=(idx > 0))

    def parseOutput(self):
        '''
        Function to parse outputs from SCRIBE.
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
            outFile = 'outFile'+str(idx)+'.csv'
            if not (workDir / outFile).exists():
                # Quit if output file does not exist
                print(str(workDir / outFile) + ' does not exist, skipping...')
                return
            OutSubDF[idx] = pd.read_csv(workDir / outFile, sep = ' ', header = None)

        # megre the dataframe by taking the maximum value from each DF
        # From here: https://stackoverflow.com/questions/20383647/pandas-selecting-by-label-sometimes-return-series-sometimes-returns-dataframe
        outDF = pd.concat(OutSubDF)
        outDF.columns= ['Gene1','Gene2','EdgeWeight']
        # Group by rows code is from here:
        # https://stackoverflow.com/questions/53114609/pandas-how-to-remove-duplicate-rows-but-keep-all-rows-with-max-value
        res = outDF[outDF['EdgeWeight'] == outDF.groupby(['Gene1','Gene2'])['EdgeWeight'].transform('max')]
        # Sort values in the dataframe
        finalDF = res.sort_values('EdgeWeight',ascending=False)

        self._write_ranked_edges(finalDF[['Gene1', 'Gene2', 'EdgeWeight']])
