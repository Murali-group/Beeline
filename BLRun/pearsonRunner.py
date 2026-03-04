import pandas as pd

from BLRun.runner import Runner


class PearsonRunner(Runner):
    """Concrete runner for pairwise Pearson correlation GRN inference.
    Runs entirely within the BEELINE conda environment; no Docker image is used.
    The image field in the config should be set to 'local'."""

    def generateInputs(self):
        '''
        Verifies that the expression data file exists in the input directory.
        No file copying is required because Pearson runs locally without Docker.

        :param self.input_dir: Path — directory containing input files
        :param self.exprData: str — expression data filename
        :raises FileNotFoundError: if the expression data file is missing
        '''
        if not (self.input_dir / self.exprData).exists():
            raise FileNotFoundError(
                f"Expression data file not found: {self.input_dir / self.exprData}")

    def run(self):
        '''
        Computes pairwise Pearson correlation between all gene pairs.
        Each gene's expression is first normalized by its maximum expression value
        across all cells, bringing values into the range [0, 1] for non-negative data.
        Genes with zero maximum expression are left unnormalized (divisor replaced with 1).
        Writes the full (genes x genes) correlation matrix to working_dir/outFile.txt.

        :param self.input_dir: Path — directory containing expression data
        :param self.exprData: str — CSV filename; rows = genes, columns = cells
        :param self.working_dir: Path — output location for outFile.txt
        :output working_dir/outFile.txt: tab-separated (genes x genes) correlation matrix
        '''
        # Read expression data: rows = genes, columns = cells
        ExpressionData = pd.read_csv(
            self.input_dir / self.exprData, header=0, index_col=0)
        if not isinstance(ExpressionData, pd.DataFrame):
            raise TypeError(f"ExpressionData must be a DataFrame, got {type(ExpressionData)}")

        # Normalize each gene (row) by its maximum expression value.
        # max=0 is replaced with 1 to avoid division by zero for silent genes.
        max_per_gene = ExpressionData.max(axis=1).replace(0, 1)
        normalized = ExpressionData.div(max_per_gene, axis=0)

        # Transpose to (cells x genes) so .corr() produces a (genes x genes) matrix.
        corr = normalized.T.corr(method='pearson')

        corr.to_csv(self.working_dir / 'outFile.txt', sep='\t')

    def parseOutput(self):
        '''
        Reads the gene x gene correlation matrix from working_dir/outFile.txt and
        writes a ranked edge list to output_dir/rankedEdges.csv.
        Both directions of each gene pair are included (the matrix is symmetric).
        Self-correlations (Gene1 == Gene2) are excluded.
        Edges are ranked by absolute correlation value, descending.

        :param self.working_dir: Path — directory containing outFile.txt
        :output output_dir/rankedEdges.csv: tab-separated edge list with columns
            Gene1 (str), Gene2 (str), EdgeWeight (float, signed Pearson r)
        '''
        outFile = self.working_dir / 'outFile.txt'
        if not outFile.exists():
            print(str(outFile) + ' does not exist, skipping...')
            return

        # Read square correlation matrix (genes x genes)
        CorrDF = pd.read_csv(outFile, sep='\t', header=0, index_col=0)
        if not isinstance(CorrDF, pd.DataFrame):
            raise TypeError(f"CorrDF must be a DataFrame, got {type(CorrDF)}")

        # Convert to long-format edge list and drop self-correlations
        stacked = CorrDF.stack().reset_index()
        stacked.columns = ['Gene1', 'Gene2', 'EdgeWeight']
        OutDF = stacked[stacked['Gene1'] != stacked['Gene2']].copy()

        # Rank by absolute correlation value, descending; retain signed EdgeWeight
        OutDF = OutDF.iloc[OutDF['EdgeWeight'].abs().argsort()[::-1]]

        self._write_ranked_edges(OutDF[['Gene1', 'Gene2', 'EdgeWeight']])
