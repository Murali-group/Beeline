from itertools import combinations
from typing import Dict, List, Set, Tuple

import pandas as pd

from BLEval.evaluator import Evaluator
from BLEval.data import EvaluationData, DatasetGroup


def _top_k_edges(ranked_edges: pd.DataFrame, k: int) -> Set[Tuple[str, str]]:
    """
    Return the set of top-k predicted edges after removing self-loops.

    Duplicate (Gene1, Gene2) pairs are deduplicated by keeping the row with
    the highest absolute EdgeWeight. A tie-aware cutoff is applied at rank k:
    bestVal = max(nonZeroMin, edgeWeightTopk) expands the selected set to
    include all predictions tied at the boundary, avoiding the large block of
    zero-weight unpredicted edges when the k-th prediction has weight 0.
    k is capped at the number of available predictions to avoid index errors.

    Parameters
    ----------
    ranked_edges : pd.DataFrame
        Predicted edge list with columns Gene1, Gene2, EdgeWeight.
    k : int
        Number of top edges to select after filtering self-loops.

    Returns
    -------
    set of (Gene1, Gene2) tuples
        The top-k predicted directed edges (expanded to include ties).
    """
    if not isinstance(ranked_edges, pd.DataFrame):
        raise TypeError(f"ranked_edges must be DataFrame, got {type(ranked_edges)}")
    if not isinstance(k, int):
        raise TypeError(f"k must be int, got {type(k)}")

    no_self = ranked_edges[ranked_edges['Gene1'] != ranked_edges['Gene2']].copy()

    if no_self.empty:
        return set()

    # Deduplicate: per (Gene1, Gene2), keep the row with the highest abs weight
    no_self['_abs'] = no_self['EdgeWeight'].abs()
    no_self = no_self.sort_values('_abs', ascending=False).drop_duplicates(
        subset=['Gene1', 'Gene2']
    )

    if no_self.empty:
        return set()

    # Cap k at the number of available predictions to avoid index errors
    maxk = min(len(no_self), k)
    edge_weight_topk = float(no_self.iloc[maxk - 1]['_abs'])

    nonzero = no_self.loc[no_self['_abs'] > 0, '_abs']
    non_zero_min = float(nonzero.min()) if not nonzero.empty else 0.0

    # Expand selection to include all ties at the boundary; avoids pulling in
    # the large block of zero-weight edges when edgeWeightTopk is 0.
    best_val = max(non_zero_min, edge_weight_topk)
    selected = no_self[no_self['_abs'] >= best_val]

    return set(zip(selected['Gene1'], selected['Gene2']))


def _jaccard(set_a: Set, set_b: Set) -> float:
    """
    Compute the Jaccard index between two sets.

    Jaccard index = |A ∩ B| / |A ∪ B|. Returns float('nan') when both sets
    are empty (the index is undefined for two empty sets).

    Parameters
    ----------
    set_a : set
        First edge set.
    set_b : set
        Second edge set.

    Returns
    -------
    float
        Jaccard index in [0, 1], or nan if both sets are empty.
    """
    if not isinstance(set_a, set):
        raise TypeError(f"set_a must be set, got {type(set_a)}")
    if not isinstance(set_b, set):
        raise TypeError(f"set_b must be set, got {type(set_b)}")

    union = set_a | set_b
    if not union:
        return float('nan')
    return len(set_a & set_b) / len(union)


def _compute_jaccard(
    run_edge_sets: Dict[str, Set[Tuple[str, str]]],
) -> Tuple[float, float]:
    """
    Compute the median and mean absolute deviation of pairwise Jaccard indices
    across a collection of run edge sets.

    All pairs of runs are compared and both the median and mean absolute
    deviation (MAD) of pairwise Jaccard scores are returned. Returns
    (nan, nan) when fewer than two runs are available (no pairs to compare).

    Parameters
    ----------
    run_edge_sets : dict[str, set]
        Mapping of run_id to the top-k predicted edge set for one algorithm.

    Returns
    -------
    tuple of (float, float)
        Median Jaccard index and mean absolute deviation across all run pairs,
        or (nan, nan) if fewer than 2 runs.
    """
    if not isinstance(run_edge_sets, dict):
        raise TypeError(f"run_edge_sets must be dict, got {type(run_edge_sets)}")

    run_ids = list(run_edge_sets.keys())

    if len(run_ids) < 2:
        return float('nan'), float('nan')

    scores: List[float] = [
        _jaccard(run_edge_sets[ri], run_edge_sets[rj])
        for ri, rj in combinations(run_ids, 2)
    ]

    s = pd.Series(scores)
    median = float(s.median())
    mad = float((s - s.mean()).abs().mean())
    return median, mad


class Jaccard(Evaluator):
    """
    Evaluator that computes the median pairwise Jaccard index of top-k
    predicted networks across runs sharing the same ground truth.

    The Jaccard index measures the stability of each algorithm's predictions:
    runs within a DatasetGroup are perturbations of the same biological system,
    so a high median Jaccard indicates that the algorithm produces consistent
    top-k edge sets regardless of sampling noise.

    k is set to the number of edges in the ground truth network (excluding
    self-loops), consistent with the EarlyPrecision evaluator. For each
    DatasetGroup, writes Jaccard.csv to dataset_path. Rows are algorithms and
    the single column is the median Jaccard index across all run pairs.
    DatasetGroups with fewer than two runs produce nan for all algorithms.
    """

    def __call__(self, evaluation_data: EvaluationData) -> None:
        """
        Compute median pairwise Jaccard per algorithm and write results to
        dataset_path/Jaccard.csv.

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
            # Determine k from the first run with an accessible ground truth.
            # All runs in a DatasetGroup share the same ground truth network,
            # so k is constant across runs.
            k: int = 0
            for run in dataset_group:
                if run.ground_truth_path.exists():
                    gt_df = self._load_ground_truth(run.ground_truth_path)
                    k = sum(
                        1 for g1, g2 in zip(gt_df['Gene1'], gt_df['Gene2'])
                        if g1 != g2
                    )
                    break

            if k == 0:
                print(
                    f"Warning: no ground truth found for dataset '{dataset_group.dataset_id}', "
                    "skipping Jaccard computation."
                )
                continue

            # Collect top-k edge sets per run per algorithm.
            # run_sets[run_id][algo] = set of top-k (Gene1, Gene2) tuples
            run_sets: Dict[str, Dict[str, Set[Tuple[str, str]]]] = {}
            for run in dataset_group:
                run_sets[run.run_id] = {
                    algo: _top_k_edges(df, k)
                    for algo, df in run.ranked_edges.items()
                }

            # Collect the union of algorithm names across all runs.
            algos = {algo for sets in run_sets.values() for algo in sets}

            if not algos:
                continue

            # For each algorithm, build a run_id → edge set mapping and compute
            # the median and MAD of pairwise Jaccard across all run pairs.
            results: Dict[str, Tuple[float, float]] = {}
            for algo in sorted(algos):
                per_run = {
                    run_id: sets[algo]
                    for run_id, sets in run_sets.items()
                    if algo in sets
                }
                results[algo] = _compute_jaccard(per_run)

            out_df = pd.DataFrame.from_dict(
                results, orient='index', columns=['MedianJaccard', 'MADJaccard']
            )
            out_df.index.name = 'Algorithm'

            dataset_group.dataset_path.mkdir(parents=True, exist_ok=True)
            out_path = dataset_group.dataset_path / 'Jaccard.csv'
            out_df.to_csv(out_path)
            print(f"Wrote Jaccard results to {out_path}")
