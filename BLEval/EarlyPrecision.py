from typing import Dict, Set, Tuple

import pandas as pd

from BLEval.evaluator import Evaluator
from BLEval.data import EvaluationData


def _compute_early_precision(
    ranked_edges: pd.DataFrame,
    true_edges: Set[Tuple[str, str]],
    k: int,
) -> float:
    """
    Compute early precision for one algorithm on one run.

    Early precision is the fraction of true positives in the top-k predicted
    edges, where k is the number of edges in the ground truth network excluding
    self-loops. Self-loops are excluded from predictions before ranking. Returns
    float('nan') if k is zero or the predicted list is empty after filtering.

    Parameters
    ----------
    ranked_edges : pd.DataFrame
        Predicted edge list with columns Gene1, Gene2, EdgeWeight.
        Higher EdgeWeight indicates greater confidence.
    true_edges : set of (str, str)
        Ground truth directed edges (self-loops may be present but do not
        affect k, which is computed by the caller before this function).
    k : int
        Number of ground truth edges excluding self-loops. Serves as both
        the cutoff rank and the denominator of the precision calculation.

    Returns
    -------
    float
        Early precision in [0, 1], or nan if the metric cannot be computed.
    """
    if not isinstance(ranked_edges, pd.DataFrame):
        raise TypeError(f"ranked_edges must be DataFrame, got {type(ranked_edges)}")
    if not isinstance(true_edges, set):
        raise TypeError(f"true_edges must be set, got {type(true_edges)}")
    if not isinstance(k, int):
        raise TypeError(f"k must be int, got {type(k)}")

    if k == 0:
        return float('nan')

    # Remove self-loops — they are never meaningful predictions
    predicted = ranked_edges[ranked_edges['Gene1'] != ranked_edges['Gene2']].copy()

    if predicted.empty:
        return float('nan')

    # Sort descending by edge weight and take the top-k predictions
    top_k = predicted.nlargest(k, 'EdgeWeight')

    true_positives = sum(
        1 for edge in zip(top_k['Gene1'], top_k['Gene2'])
        if edge in true_edges
    )

    return true_positives / k


class EarlyPrecision(Evaluator):
    """
    Evaluator that computes early precision for each algorithm against the
    ground truth network.

    Early precision is defined as the fraction of true positives among the
    top-k predicted edges, where k equals the number of edges in the ground
    truth network (excluding self-loops). For each DatasetGroup, writes
    EarlyPrecision.csv to dataset_path. Rows are algorithms and columns are
    run_ids. Runs whose ground truth file is missing are skipped with a warning.
    """

    def __call__(self, evaluation_data: EvaluationData) -> None:
        """
        Compute early precision per algorithm per run and write results to
        dataset_path/EarlyPrecision.csv.

        Parameters
        ----------
        evaluation_data : EvaluationData
            Loaded predicted networks organised by dataset and run.

        Returns
        -------
        None
        """
        if not isinstance(evaluation_data, EvaluationData):
            raise TypeError(f"evaluation_data must be EvaluationData, got {type(evaluation_data)}")

        for dataset_group in evaluation_data:
            # results[algo][run_id] = early precision score
            results: Dict[str, Dict[str, float]] = {}

            for run in dataset_group:
                if not run.ground_truth_path.exists():
                    print(
                        f"Warning: ground truth not found at {run.ground_truth_path}, "
                        f"skipping run '{run.run_id}'."
                    )
                    continue

                true_edges = self._load_ground_truth(run.ground_truth_path)

                # k = number of ground truth edges excluding self-loops
                k = sum(1 for g1, g2 in true_edges if g1 != g2)

                for algo, ranked_edges_df in run.ranked_edges.items():
                    score = _compute_early_precision(ranked_edges_df, true_edges, k)
                    results.setdefault(algo, {})[run.run_id] = score

            if not results:
                continue

            # Build output DataFrame: rows = algorithms, columns = run_ids
            out_df = pd.DataFrame(results).T
            out_df.index.name = 'Algorithm'

            dataset_group.dataset_path.mkdir(parents=True, exist_ok=True)
            out_path = dataset_group.dataset_path / 'EarlyPrecision.csv'
            out_df.to_csv(out_path)
            print(f"Wrote EarlyPrecision results to {out_path}")
