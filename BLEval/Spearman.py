from itertools import combinations
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from BLEval.evaluator import Evaluator
from BLEval.data import EvaluationData, RunResult


def _align_weights(
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
) -> tuple:
    """
    Align two ranked-edge DataFrames on their union of non-self-loop edges.

    Edges present in one DataFrame but absent from the other receive a weight
    of 0 (not predicted = lowest confidence). Returns two equal-length numpy
    arrays of EdgeWeights in the same edge order, suitable for rank correlation.

    Parameters
    ----------
    df_a : pd.DataFrame
        Ranked edges from run A with columns Gene1, Gene2, EdgeWeight.
    df_b : pd.DataFrame
        Ranked edges from run B with columns Gene1, Gene2, EdgeWeight.

    Returns
    -------
    weights_a : np.ndarray
        EdgeWeight values for each edge in the union, from run A.
    weights_b : np.ndarray
        EdgeWeight values for each edge in the union, from run B.
    """
    if not isinstance(df_a, pd.DataFrame):
        raise TypeError(f"df_a must be DataFrame, got {type(df_a)}")
    if not isinstance(df_b, pd.DataFrame):
        raise TypeError(f"df_b must be DataFrame, got {type(df_b)}")

    def _to_series(df: pd.DataFrame) -> pd.Series:
        """Return a Series of EdgeWeight indexed by (Gene1, Gene2), self-loops removed.

        Duplicate (Gene1, Gene2) pairs are deduplicated by keeping the row with
        the highest EdgeWeight, matching the top-k selection used elsewhere.
        """
        no_loops = df[df['Gene1'] != df['Gene2']].copy()
        # Keep highest-weight entry per directed edge to ensure a unique index.
        no_loops = no_loops.sort_values('EdgeWeight', ascending=False)
        no_loops = no_loops.drop_duplicates(subset=['Gene1', 'Gene2'], keep='first')
        idx = pd.MultiIndex.from_arrays([no_loops['Gene1'], no_loops['Gene2']])
        return pd.Series(no_loops['EdgeWeight'].values, index=idx)

    s_a = _to_series(df_a)
    s_b = _to_series(df_b)

    # Reindex both to the union; missing edges default to 0
    union_idx = s_a.index.union(s_b.index)
    weights_a = s_a.reindex(union_idx, fill_value=0.0).values
    weights_b = s_b.reindex(union_idx, fill_value=0.0).values

    return weights_a, weights_b


def _compute_median_spearman(
    runs: List[RunResult],
    algo: str,
) -> float:
    """
    Compute the median pairwise Spearman correlation for one algorithm across runs.

    All pairs of runs that both have predictions for the given algorithm are
    evaluated. The Spearman correlation is computed on the union of edges
    (missing edges filled with weight 0). Returns float('nan') when fewer than
    two runs have predictions for this algorithm, or all pairwise correlations
    are undefined.

    Parameters
    ----------
    runs : list of RunResult
        All runs within a DatasetGroup.
    algo : str
        Algorithm name to evaluate.

    Returns
    -------
    float
        Median Spearman correlation across all pairwise run combinations,
        or nan if the metric cannot be computed.
    """
    if not isinstance(runs, list):
        raise TypeError(f"runs must be list, got {type(runs)}")
    if not isinstance(algo, str):
        raise TypeError(f"algo must be str, got {type(algo)}")

    # Only keep runs that have predictions for this algorithm
    algo_runs = [run for run in runs if algo in run.ranked_edges]

    if len(algo_runs) < 2:
        return float('nan')

    correlations: List[float] = []
    for run_a, run_b in combinations(algo_runs, 2):
        weights_a, weights_b = _align_weights(
            run_a.ranked_edges[algo],
            run_b.ranked_edges[algo],
        )

        # Spearman is undefined when all weights are identical in either vector
        if np.all(weights_a == weights_a[0]) or np.all(weights_b == weights_b[0]):
            continue

        corr, _ = spearmanr(weights_a, weights_b)
        if not np.isnan(corr):
            correlations.append(float(corr))

    if not correlations:
        return float('nan')

    return float(np.median(correlations))


class Spearman(Evaluator):
    """
    Evaluator that measures how consistently each algorithm ranks edges across
    runs within a DatasetGroup.

    For each algorithm, all pairs of runs are compared by computing the Spearman
    rank correlation of their EdgeWeight vectors (edges absent from a run receive
    weight 0). The median of all pairwise correlations is reported. For each
    DatasetGroup, writes Spearman.csv to dataset_path. Rows are algorithms and
    the single column MedianSpearman holds the median correlation. Algorithms
    appearing in fewer than two runs receive nan.
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

            # Collect all algorithm names that appear in at least one run
            algos = sorted({algo for run in runs for algo in run.ranked_edges})

            if not algos:
                continue

            # results[algo] = median Spearman correlation
            results: Dict[str, float] = {}
            for algo in algos:
                results[algo] = _compute_median_spearman(runs, algo)

            out_df = pd.DataFrame.from_dict(
                results, orient='index', columns=['MedianSpearman']
            )
            out_df.index.name = 'Algorithm'

            dataset_group.dataset_path.mkdir(parents=True, exist_ok=True)
            out_path = dataset_group.dataset_path / 'Spearman.csv'
            out_df.to_csv(out_path)
            print(f"Wrote Spearman results to {out_path}")
