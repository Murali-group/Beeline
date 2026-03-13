from itertools import permutations
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from BLEval.evaluator import Evaluator
from BLEval.data import EvaluationData, RunResult


def _to_series(df: pd.DataFrame) -> pd.Series:
    """
    Convert a ranked-edge DataFrame to a signed-weight Series indexed by (Gene1, Gene2).

    Self-loops are removed. Duplicate (Gene1, Gene2) pairs are deduplicated by
    keeping the row with the highest absolute EdgeWeight, preserving its sign.

    Parameters
    ----------
    df : pd.DataFrame
        Ranked edges with columns Gene1, Gene2, EdgeWeight.

    Returns
    -------
    pd.Series
        Signed EdgeWeight values indexed by (Gene1, Gene2) MultiIndex.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"df must be DataFrame, got {type(df)}")

    no_loops = df[df['Gene1'] != df['Gene2']].copy()
    no_loops['_abs'] = no_loops['EdgeWeight'].abs()
    no_loops = no_loops.sort_values('_abs', ascending=False).drop(columns=['_abs'])
    no_loops = no_loops.drop_duplicates(subset=['Gene1', 'Gene2'], keep='first')
    idx = pd.MultiIndex.from_arrays([no_loops['Gene1'], no_loops['Gene2']])
    return pd.Series(no_loops['EdgeWeight'].values, index=idx)


def _compute_median_spearman(
    runs: List[RunResult],
    algo: str,
    gene_set: set,
) -> Tuple[float, float]:
    """
    Compute the median and mean absolute deviation of pairwise Spearman
    correlations for one algorithm across runs.

    Weight arrays for each run are precomputed once against the fixed edge
    universe, then all pairwise correlations are obtained in a single
    scipy.stats.spearmanr matrix call rather than one call per pair. Runs
    whose weight vector is constant (all weights identical) are excluded before
    the matrix call, since Spearman is undefined for constant inputs. Returns
    (nan, nan) when fewer than two usable runs exist or all correlations are
    undefined.

    Parameters
    ----------
    runs : list of RunResult
        All runs within a DatasetGroup.
    algo : str
        Algorithm name to evaluate.
    gene_set : set
        Unique genes from the ground truth network. Determines the fixed edge
        universe (all directed non-self-loop pairs) used for alignment.

    Returns
    -------
    tuple of (float, float)
        Median Spearman correlation and mean absolute deviation across all
        pairwise run combinations, or (nan, nan) if the metric cannot be computed.
    """
    if not isinstance(runs, list):
        raise TypeError(f"runs must be list, got {type(runs)}")
    if not isinstance(algo, str):
        raise TypeError(f"algo must be str, got {type(algo)}")
    if not isinstance(gene_set, set):
        raise TypeError(f"gene_set must be set, got {type(gene_set)}")

    algo_runs = [run for run in runs if algo in run.ranked_edges]

    if len(algo_runs) < 2:
        return float('nan'), float('nan')

    # Build fixed edge universe once — all directed pairs among GT genes
    possible_edges = list(permutations(sorted(gene_set), 2))
    if not possible_edges:
        return float('nan'), float('nan')
    fixed_idx = pd.MultiIndex.from_tuples(possible_edges)

    # Precompute weight array for each run once, skipping constant-weight runs
    # (Spearman is undefined when one vector has zero variance).
    weight_arrays: List[np.ndarray] = []
    for run in algo_runs:
        w = _to_series(run.ranked_edges[algo]).reindex(fixed_idx, fill_value=0.0).values
        if not np.all(w == w[0]):
            weight_arrays.append(w)

    if len(weight_arrays) < 2:
        return float('nan'), float('nan')

    # Stack into (n_edges, n_runs) matrix and compute all pairwise correlations
    # in one call; much faster than calling spearmanr once per pair.
    matrix = np.column_stack(weight_arrays)
    result = spearmanr(matrix)

    # spearmanr returns a scalar correlation when there are exactly 2 columns,
    # and a correlation matrix when there are more.
    # Use .correlation (scipy < 1.7) which is available in all supported versions.
    n = len(weight_arrays)
    if n == 2:
        corr_arr = np.array([float(result.correlation)])
    else:
        corr_matrix = np.asarray(result.correlation)
        rows, cols = np.triu_indices(n, k=1)
        corr_arr = corr_matrix[rows, cols]

    corr_arr = corr_arr[~np.isnan(corr_arr)]

    if corr_arr.size == 0:
        return float('nan'), float('nan')

    median = float(np.median(corr_arr))
    mad = float(np.mean(np.abs(corr_arr - np.mean(corr_arr))))
    return median, mad


class Spearman(Evaluator):
    """
    Evaluator that measures how consistently each algorithm ranks edges across
    runs within a DatasetGroup.

    For each algorithm, all pairs of runs are compared by computing the Spearman
    rank correlation of their signed EdgeWeight vectors over the fixed universe
    of all directed gene pairs in the ground truth (missing edges receive weight
    0, signed weights preserved). Weight arrays are precomputed once per run
    and all pairwise correlations are obtained via a single matrix call to
    scipy.stats.spearmanr. The median and mean absolute deviation of all
    pairwise correlations are reported. For each DatasetGroup, writes
    Spearman.csv to dataset_path. Rows are algorithms; columns are
    MedianSpearman and MADSpearman.
    """

    def __call__(self, evaluation_data: EvaluationData) -> None:
        """
        Compute median pairwise Spearman correlation per algorithm and write
        results to dataset_path/Spearman.csv.

        Parameters
        ----------
        evaluation_data : EvaluationData
            Loaded predicted networks organised by dataset and run.

        Returns
        -------
        None
        """
        if not isinstance(evaluation_data, EvaluationData):
            raise TypeError(
                f"evaluation_data must be EvaluationData, got {type(evaluation_data)}"
            )

        for dataset_group in evaluation_data:
            runs = dataset_group.runs

            algos = sorted({algo for run in runs for algo in run.ranked_edges})

            if not algos:
                continue

            # Load ground truth gene set from the first run with an existing file.
            # All runs in a DatasetGroup share the same ground truth network.
            gt_path = next(
                (run.ground_truth_path for run in runs
                 if run.ground_truth_path.exists()),
                None,
            )
            if gt_path is None:
                print(
                    f"Warning: no ground truth found for dataset "
                    f"'{dataset_group.dataset_id}', skipping Spearman."
                )
                continue

            gt_df = pd.read_csv(gt_path, header=0)
            gene_set = set(gt_df['Gene1']).union(set(gt_df['Gene2']))

            # results[algo] = (median, mad)
            results: Dict[str, Tuple[float, float]] = {}
            for algo in algos:
                results[algo] = _compute_median_spearman(runs, algo, gene_set)

            out_df = pd.DataFrame.from_dict(
                results, orient='index', columns=['MedianSpearman', 'MADSpearman']
            )
            out_df.index.name = 'Algorithm'

            dataset_group.dataset_path.mkdir(parents=True, exist_ok=True)
            out_path = dataset_group.dataset_path / 'Spearman.csv'
            out_df.to_csv(out_path)
            print(f"Wrote Spearman results to {out_path}")
