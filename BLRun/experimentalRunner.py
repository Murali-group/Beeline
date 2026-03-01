import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import wasserstein_distance

from BLRun.runner import Runner

class ExperimentalRunner(Runner):
    """Concrete runner for experimental GRN inference algorithms.
    Runs entirely within the BEELINE conda environment; no Docker image is used.
    The image field in the config should be set to 'local'."""

    def __init__(self, root: Path, config: dict):
        super().__init__(root, config)

        # Number of equal-width bins used to discretise the [0, 1] pseudotime axis
        # when building per-gene expression distributions.
        self._N_BINS = 100

        # Scale parameter (σ_dir) for the sigmoid directional-bias term in parseOutput.
        # Controls how quickly the directional preference saturates as the centre-of-mass
        # shift between two genes grows.  Expressed as a fraction of the [0, 1]
        # pseudotime range; 0.1 means a 10 % shift gives sigmoid(1) ≈ 0.73 weight.
        self._LAG_SCALE = 0.1

        # Scale parameter (σ_prox) for the Gaussian proximity-bias term in parseOutput.
        # Controls how quickly scores are suppressed as the centre-of-mass distance grows.
        # Expressed as a fraction of the [0, 1] pseudotime range; 0.3 means a 30 % shift
        # gives exp(-0.5) ≈ 0.61 weight, and a 60 % shift gives exp(-2) ≈ 0.14 weight.
        self._PROXIMITY_SCALE = 0.1

    def generateInputs(self):
        '''
        Verifies that the expression data and pseudotime files exist in the
        input directory. No file copying is required because this runner runs
        locally without Docker.

        :param self.input_dir: Path — directory containing input files
        :param self.exprData: str — expression data filename
        :param self.pseudoTimeData: str — pseudotime data filename
        :raises FileNotFoundError: if either input file is missing
        '''
        for fname in (self.exprData, self.pseudoTimeData):
            if not (self.input_dir / fname).exists():
                raise FileNotFoundError(
                    f"Input file not found: {self.input_dir / fname}")

    def run(self):
        '''
        Per-gene pseudotime expression distributions (working_dir/distributions.txt):
           For each gene, cells are binned along the pseudotime axis (_N_BINS equal-width
           bins over [0, 1]). Expression is normalised per gene to the fraction of total
           expression across all cells before accumulation; the resulting histogram is
           further scaled so its area (sum × bin_width) equals 1.0.
           Rows = genes, columns = bin-centre pseudotime values.

        :param self.input_dir: Path — directory containing expression and pseudotime data
        :param self.exprData: str — CSV filename; rows = genes, columns = cells
        :param self.pseudoTimeData: str — CSV filename; rows = cells, column = PseudoTime
        :param self.working_dir: Path — output location
        :output working_dir/distributions.txt: tab-separated (genes x bins) distribution matrix
        '''
        # -- Load data --------------------------------------------------------
        ExpressionData = pd.read_csv(
            self.input_dir / self.exprData, header=0, index_col=0)
        if not isinstance(ExpressionData, pd.DataFrame):
            raise TypeError(f"ExpressionData must be a DataFrame, got {type(ExpressionData)}")

        # PseudoTime: rows = cells, single column
        PseudoTime = pd.read_csv(
            self.input_dir / self.pseudoTimeData, index_col=0).iloc[:, 0]
        if not isinstance(PseudoTime, pd.Series):
            raise TypeError(f"PseudoTime must be a Series, got {type(PseudoTime)}")

        # Align to cells present in both datasets
        common_cells = ExpressionData.columns.intersection(PseudoTime.index)
        ExpressionData = ExpressionData[common_cells]
        PseudoTime     = PseudoTime[common_cells]

        # -- Pseudotime normalisation to [0, 1] -------------------------------
        pt_min, pt_max = PseudoTime.min(), PseudoTime.max()
        if pt_max - pt_min > 0:
            pt_norm = (PseudoTime - pt_min) / (pt_max - pt_min)
        else:
            # Degenerate case: all cells share the same pseudotime value.
            pt_norm = pd.Series(0.0, index=PseudoTime.index)

        # -- Per-gene expression normalisation to fraction of total -----------
        # gene_totals[g] = total expression of gene g summed over all cells.
        # Silent genes (total = 0) are left with a divisor of 1 to avoid NaN.
        gene_totals = ExpressionData.sum(axis=1).replace(0, 1)
        expr_frac   = ExpressionData.div(gene_totals, axis=0)

        # -- Pseudotime distribution for each gene ----------------------------
        bin_edges  = np.linspace(0, 1, self._N_BINS + 1)
        bin_width  = 1.0 / self._N_BINS
        bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2

        pt_values = pt_norm.values  # (n_cells,) array, aligned with expr_frac columns

        dist_rows = {}
        for gene in ExpressionData.index:
            weights = expr_frac.loc[gene].values
            hist, _ = np.histogram(pt_values, bins=bin_edges, weights=weights)
            # Normalise so total area (sum × bin_width) = 1.0.
            area = hist.sum() * bin_width
            if area > 0:
                hist = hist / area
            dist_rows[gene] = hist

        dist_df = pd.DataFrame.from_dict(dist_rows, orient='index', columns=bin_centres)
        dist_df.to_csv(self.working_dir / 'distributions.txt', sep='\t')

    def parseOutput(self):
        '''
        Reads the per-gene pseudotime distributions from working_dir/distributions.txt
        and writes an asymmetric ranked edge list to output_dir/rankedEdges.csv.

        Edge weight formula for directed edge A → B:
            w(A→B) = sign(r) × |r| × sigmoid((μ_B − μ_A) / σ_dir)
                                    × exp(−(μ_B − μ_A)² / 2σ_prox²)

        where r is the Pearson correlation of the two density vectors and μ_g is the
        centre-of-mass (weighted-mean pseudotime) of gene g's distribution.

        The sigmoid term (σ_dir = _LAG_SCALE) is directional: when A's centre of mass
        precedes B's, it is > 0.5, boosting A→B relative to B→A.

        The Gaussian term (σ_prox = _PROXIMITY_SCALE) is symmetric: it boosts pairs
        whose centres of mass are close in pseudotime and suppresses pairs that are
        far apart, regardless of direction.

        Self-edges (Gene1 == Gene2) are excluded.
        Edges are ranked by absolute EdgeWeight, descending.

        :param self.working_dir: Path — directory containing distributions.txt
        :output output_dir/rankedEdges.csv: tab-separated edge list with columns
            Gene1 (str), Gene2 (str), EdgeWeight (float, asymmetric lag-biased score)
        '''
        distFile = self.working_dir / 'distributions.txt'
        if not distFile.exists():
            print(str(distFile) + ' does not exist, skipping...')
            return

        dist_df = pd.read_csv(distFile, sep='\t', header=0, index_col=0)
        if not isinstance(dist_df, pd.DataFrame):
            raise TypeError(f"dist_df must be a DataFrame, got {type(dist_df)}")

        bin_centres = dist_df.columns.values.astype(float)  # (n_bins,)
        bin_width   = bin_centres[1] - bin_centres[0]

        # Centre of mass: expected pseudotime under each gene's density.
        # mu[i] = Σ(bin_centre × density × bin_width) for gene i.
        mu = (dist_df.values * bin_centres[np.newaxis, :] * bin_width).sum(axis=1)

        # Symmetric Pearson correlation matrix (genes × genes) from density vectors.
        r_matrix = dist_df.T.corr(method='pearson').values

        # Shift matrix: shift[i, j] = μ_j − μ_i.
        # Positive shift means gene i's expression precedes gene j's → supports i→j.
        shift_matrix = mu[np.newaxis, :] - mu[:, np.newaxis]

        # Directional bias: sigmoid > 0.5 when shift > 0 (i leads j), < 0.5 otherwise.
        lag_bias = 1.0 / (1.0 + np.exp(-shift_matrix / self._LAG_SCALE))

        # Proximity bias: Gaussian on |shift|, symmetric between both directions.
        # Boosts pairs with close centres of mass; suppresses distant pairs.
        proximity_bias = np.exp(-shift_matrix ** 2 / (2.0 * self._PROXIMITY_SCALE ** 2))

        # Asymmetric edge weight: sign(r) × |r| × directional_bias × proximity_bias.
        weighted = np.sign(r_matrix) * np.abs(r_matrix) * lag_bias * proximity_bias

        genes       = dist_df.index.tolist()
        weighted_df = pd.DataFrame(weighted, index=genes, columns=genes)

        # Convert to long-format edge list, drop self-edges, rank by |weight|.
        stacked = weighted_df.stack().reset_index()
        stacked.columns = ['Gene1', 'Gene2', 'EdgeWeight']
        OutDF = stacked[stacked['Gene1'] != stacked['Gene2']].copy()
        OutDF = OutDF.reindex(OutDF['EdgeWeight'].abs().sort_values(ascending=False).index)

        self._write_ranked_edges(OutDF[['Gene1', 'Gene2', 'EdgeWeight']])


class GeneticRunner(Runner):
    """Concrete runner for the GENETIC GRN inference algorithm.
    Runs entirely within the BEELINE conda environment; no Docker image is used.
    The image field in the config should be set to 'local'."""

    def generateInputs(self):
        '''
        Verifies that the expression data file exists in the input directory.
        No file copying is required because GENETIC runs locally without Docker.

        :param self.input_dir: Path — directory containing input files
        :param self.exprData: str — expression data filename
        :raises FileNotFoundError: if the expression data file is missing
        '''
        if not (self.input_dir / self.exprData).exists():
            raise FileNotFoundError(
                f"Expression data file not found: {self.input_dir / self.exprData}")

    def run(self):
        '''
        Constructs an equal-weight fully-connected DAG over all genes, simulates
        gene expression via a nonlinear additive noise model (ANM), computes
        per-gene Wasserstein loss, then attributes that loss to individual edges
        via causal intervention (edge ablation).

        DAG construction: the upper triangle of the adjacency matrix is used
        (adj[i, j] = 1.0 for i < j, 0.0 otherwise), so gene i can only regulate
        gene j when i precedes j in index order. This guarantees acyclicity and
        makes topological order identical to gene index order.

        Simulation (nonlinear SEM):
            - Root genes (no parents): x_i = ε_i,  ε_i ~ N(0, 1)
            - Non-root genes: x_i = sigmoid(Σ_j w_ji · x_j) + ε_i,  ε_i ~ N(0, noise_scale²)
        Noise is pre-generated once and reused for all ablation simulations so
        that attribution differences reflect edge removal only.

        Edge attribution and weight update:
            For each edge i → j, ablate it and re-simulate genes j..n-1 with
            the same fixed noise. attribution[i, j] is the change in Wasserstein
            distance for gene j. Positive = edge was beneficial; negative = harmful.
            Weights are updated without clipping, so negative weights are allowed
            and represent inhibitory regulatory interactions.

        :param self.input_dir: Path — directory containing expression data
        :param self.exprData: str — CSV filename; rows = genes, columns = cells
        :param self.working_dir: Path — output location
        :param self.params['n_cells']: str (optional) — number of cells to simulate; default 100
        :param self.params['noise_scale']: str (optional) — std dev of additive noise; default 0.1
        :param self.params['lambda']: str (optional) — step size for edge weight update; default 0.1
        :param self.params['n_iterations']: str (optional) — number of update iterations; default 100
        :output working_dir/outFile.txt: tab-separated (genes x genes) updated DAG adjacency matrix
        :output working_dir/simulated_expression.txt: tab-separated (genes x cells) expression matrix
        :output working_dir/loss.txt: tab-separated (gene, WassersteinDistance) per-gene loss
        :output working_dir/edge_attribution.txt: tab-separated (genes x genes) attribution matrix
        '''
        # Read expression data to obtain the gene list
        ExpressionData = pd.read_csv(
            self.input_dir / self.exprData, header=0, index_col=0)
        if not isinstance(ExpressionData, pd.DataFrame):
            raise TypeError(f"ExpressionData must be a DataFrame, got {type(ExpressionData)}")

        genes = ExpressionData.index.tolist()
        n = len(genes)

        # -- DAG adjacency matrix (upper triangle of fully-connected graph) ---
        # adj[i, j] = 1.0 iff i < j; self-edges and back-edges are 0.
        adj = np.triu(np.ones((n, n), dtype=float), k=1)
        adj_df = pd.DataFrame(adj, index=genes, columns=genes)
        adj_df.to_csv(self.working_dir / 'outFile.txt', sep='\t')

        # -- Simulation parameters (always strings in self.params) ------------
        n_cells      = int(self.params.get('n_cells', 100))
        noise_scale  = float(self.params.get('noise_scale', 0.1))
        lam          = float(self.params.get('lambda', 0.1))
        n_iterations = int(self.params.get('n_iterations', 100))

        real = ExpressionData.values  # (n_genes, n_real_cells); fixed across all iterations

        # -- Iterative simulate → loss → attribute → update -------------------
        # Each iteration draws fresh noise (diverse cell sampling across iterations)
        # but holds it fixed within the iteration so ablation comparisons are fair.
        # adj evolves in-place; expr/per_gene_distances/attribution from the final
        # iteration are written to disk after the loop.
        expr               = None
        per_gene_distances = None
        attribution        = None

        for _iter in range(n_iterations):

            # Fresh noise each iteration; reused within the iteration.
            # noise : ndarray (n_genes, n_cells)
            noise = np.random.randn(n, n_cells)

            # -- Nonlinear SEM simulation -------------------------------------
            # expr[i, k] = expression of gene i in cell k.
            # Topological order = index order because adj is upper-triangular.
            expr = np.zeros((n, n_cells))
            for i in range(n):
                parent_mask = adj[:, i] != 0  # parents of gene i
                if not parent_mask.any():
                    expr[i, :] = noise[i, :]
                else:
                    signal    = adj[parent_mask, i] @ expr[parent_mask, :]
                    activated = 1.0 / (1.0 + np.exp(-signal))
                    expr[i, :] = activated + noise[i, :] * noise_scale

            # -- Wasserstein loss per gene ------------------------------------
            per_gene_distances = np.array([
                wasserstein_distance(real[i, :], expr[i, :])
                for i in range(n)
            ])
            if not isinstance(per_gene_distances, np.ndarray):
                raise TypeError(
                    f"per_gene_distances must be an ndarray, got {type(per_gene_distances)}")

            # -- Edge attribution via intervention ----------------------------
            # attribution[i, j] = W(real[j], sim_ablated[j]) − W(real[j], sim_original[j])
            # Positive: edge beneficial. Negative: edge harmful.
            # attribution : ndarray (n_genes, n_genes)
            attribution = np.zeros((n, n), dtype=float)
            for j in range(n):
                for i in np.where(adj[:, j] != 0)[0]:
                    adj_ablated       = adj.copy()
                    adj_ablated[i, j] = 0.0
                    expr_ablated      = expr.copy()
                    for k in range(j, n):
                        pm = adj_ablated[:, k] != 0
                        if not pm.any():
                            expr_ablated[k, :] = noise[k, :]
                        else:
                            sig = adj_ablated[pm, k] @ expr_ablated[pm, :]
                            expr_ablated[k, :] = 1.0 / (1.0 + np.exp(-sig)) + noise[k, :] * noise_scale
                    attribution[i, j] = wasserstein_distance(real[j, :], expr_ablated[j, :]) - per_gene_distances[j]

            # -- Edge weight update -------------------------------------------
            # adj[i, j] += lam * attribution[i, j].
            # No clipping: negative weights represent inhibitory interactions.
            adj = adj + lam * attribution

        # -- Write final-iteration outputs ------------------------------------
        cell_names = [f'Cell_{k}' for k in range(n_cells)]

        expr_df = pd.DataFrame(expr, index=genes, columns=cell_names)
        if not isinstance(expr_df, pd.DataFrame):
            raise TypeError(f"expr_df must be a DataFrame, got {type(expr_df)}")
        expr_df.to_csv(self.working_dir / 'simulated_expression.txt', sep='\t')

        loss_df = pd.Series(per_gene_distances, index=genes, name='WassersteinDistance')
        if not isinstance(loss_df, pd.Series):
            raise TypeError(f"loss_df must be a Series, got {type(loss_df)}")
        loss_df.to_csv(self.working_dir / 'loss.txt', sep='\t', header=True)

        attribution_df = pd.DataFrame(attribution, index=genes, columns=genes)
        if not isinstance(attribution_df, pd.DataFrame):
            raise TypeError(
                f"attribution_df must be a DataFrame, got {type(attribution_df)}")
        attribution_df.to_csv(self.working_dir / 'edge_attribution.txt', sep='\t')

        adj_df = pd.DataFrame(adj, index=genes, columns=genes)
        if not isinstance(adj_df, pd.DataFrame):
            raise TypeError(f"adj_df must be a DataFrame, got {type(adj_df)}")
        adj_df.to_csv(self.working_dir / 'outFile.txt', sep='\t')

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
            print(str(outFile) + 'does not exist, skipping...')
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
