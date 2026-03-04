import os
import pandas as pd
import numpy as np

from BLRun.runner import Runner


class SCODERunner(Runner):
    """Concrete runner for the SCODE GRN inference algorithm."""

    def generateInputs(self):
        '''
        Function to generate desired inputs for SCODE.
        If the folder/files under self.input_dir exist,
        this function will not do anything.
        '''

        ExpressionData = pd.read_csv(self.input_dir / self.exprData,
                                         header = 0, index_col = 0)
        PTData = pd.read_csv(self.input_dir / self.pseudoTimeData,
                             header = 0, index_col = 0)

        colNames = PTData.columns
        for idx in range(len(colNames)):
            # Create output subdirectory in advance to prevent docker from
            # creating it with root-exclusive permissions
            (self.working_dir / str(idx)).mkdir(exist_ok=True)

            # Select cells belonging to each pseudotime trajectory
            colName = colNames[idx]
            index = PTData[colName].index[PTData[colName].notnull()]
            exprName = "ExpressionData"+str(idx)+".csv"
            ExpressionData.loc[:,index].to_csv(self.working_dir / exprName,
                                     sep = '\t', header  = False, index = False)
            cellName = "PseudoTime"+str(idx)+".csv"
            ptDF = PTData.loc[index,[colName]]
            # SCODE expects a column labeled PseudoTime.
            ptDF.rename(columns = {colName:'PseudoTime'}, inplace = True)
            # output file
            ptDF.to_csv(self.working_dir / cellName,
                                     sep = '\t', header  = False)

    def run(self):
        '''
        Function to run SCODE algorithm
        '''

        z = str(self.params['z'])
        nIter = str(self.params['nIter'])
        nRep = str(self.params['nRep'])

        PTData = pd.read_csv(self.input_dir / self.pseudoTimeData,
                             header = 0, index_col = 0)

        colNames = PTData.columns

        for idx in range(len(colNames)):

            ExpressionData = pd.read_csv(self.working_dir /
                                         ("ExpressionData"+str(idx)+".csv"),
                                         header = None, index_col = None, sep ='\t')
            nCells = str(ExpressionData.shape[1])
            nGenes = str(ExpressionData.shape[0])

            cmdToRun = ' '.join(['docker run --rm',
                                f'--user {os.getuid()}:{os.getgid()}',
                                '-e HOME=/tmp',
                                f"-v {self.working_dir}:/usr/working_dir",
                                f'{self.image} /bin/sh -c \"time -v -o',
                                "/usr/working_dir/time" + str(idx) + ".txt",
                                'ruby run_R.rb',
                                "/usr/working_dir/ExpressionData" + str(idx) + ".csv",
                                "/usr/working_dir/PseudoTime" + str(idx) + ".csv",
                                "/usr/working_dir/" + str(idx),
                                nGenes, z, nCells, nIter, nRep, '\"'])

            self._run_docker(cmdToRun, append=(idx > 0))

    def parseOutput(self):
        '''
        Function to parse outputs from SCODE.
        '''
        workDir = self.working_dir

        PTData = pd.read_csv(self.input_dir / self.pseudoTimeData,
                             header = 0, index_col = 0)
        colNames = PTData.columns
        for indx in range(len(colNames)):
            # Read output
            outFile = str(indx)+'/meanA.txt'
            if not (workDir / outFile).exists():
                # Quit if output file does not exist
                print(str(workDir / outFile) + ' does not exist, skipping...')
                return
            OutDF = pd.read_csv(workDir / outFile, sep = '\t', header = None)

            # Sort values in a matrix using code from:
            # https://stackoverflow.com/questions/21922806/sort-values-of-matrix-in-python
            OutMatrix = np.abs(OutDF.values)
            idx = np.argsort(OutMatrix, axis = None)[::-1]
            rows, cols = np.unravel_index(idx, OutDF.shape)
            DFSorted = OutMatrix[rows, cols]

            # read input file for list of gene names
            ExpressionData = pd.read_csv(self.input_dir / self.exprData,
                                             header = 0, index_col = 0)
            GeneList = list(ExpressionData.index)

            outFile = open(workDir / ('outFile'+str(indx)+'.csv'),'w')
            outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

            for row, col, val in zip(rows, cols, DFSorted):
                outFile.write('\t'.join([GeneList[row],GeneList[col],str(val)])+'\n')
            outFile.close()

        OutSubDF = [0]*len(colNames)
        for indx in range(len(colNames)):
            outFile = 'outFile'+str(indx)+'.csv'
            OutSubDF[indx] = pd.read_csv(workDir / outFile, sep = '\t', header = 0)

            OutSubDF[indx].EdgeWeight = np.abs(OutSubDF[indx].EdgeWeight)

        outDF = pd.concat(OutSubDF)
        FinalDF = outDF[outDF['EdgeWeight'] == outDF.groupby(['Gene1','Gene2'])['EdgeWeight'].transform('max')]
        FinalDF = FinalDF.sort_values(['EdgeWeight'], ascending=False)
        self._write_ranked_edges(FinalDF[['Gene1', 'Gene2', 'EdgeWeight']])
