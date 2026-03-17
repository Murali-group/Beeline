from itertools import permutations, product
from typing import Dict, Set, Tuple

import pandas as pd

from BLEval.evaluator import Evaluator
from BLEval.data import EvaluationData


def _build_edge_universe(
    gt_df: pd.DataFrame,
    tf_edges: bool,
) -> Tuple[Set[Tuple[str, str]], Set[Tuple[str, str]]]:
    """
    Build the candidate edge universe and the set of true edges from a GT DataFrame.

    When tf_edges is False the universe is all directed non-self-loop pairs
    among GT genes. When tf_edges is True the universe is restricted to edges
    whose source gene appears as Gene1 in the ground truth (i.e. TF→gene edges),
    which is appropriate for experimental scRNA-seq data where regulators are
    known. True edges are the intersection of GT edges with the universe so
    that the random-baseline computation is consistent with the candidate pool.

    Parameters
    ----------
    gt_df : pd.DataFrame
        Ground truth network with at least columns Gene1, Gene2. Self-loops
        and duplicate rows are removed internally.
    tf_edges : bool
        When True, restrict the universe to TF→gene edges (no self-loops).
        When False, use all directed non-self-loop gene pairs.

    Returns
    -------
    possible_edges : set of (str, str)
        Complete candidate edge universe (no self-loops).
    true_edges : set of (str, str)
        Ground truth edges within the candidate universe (no self-loops).
    """
    if not isinstance(gt_df, pd.DataFrame):
        raise TypeError(f"gt_df must be DataFrame, got {type(gt_df)}")
    if not isinstance(tf_edges, bool):
        raise TypeError(f"tf_edges must be bool, got {type(tf_edges)}")

    gt_no_self = gt_df[gt_df['Gene1'] != gt_df['Gene2']].drop_duplicates()
    unique_nodes = sorted(set(gt_df['Gene1']).union(set(gt_df['Gene2'])))

    if tf_edges:
        tf_genes = set(gt_no_self['Gene1'])
        all_genes = set(unique_nodes)
        # All TF→gene edges: source is a known TF, target is any gene, no self-loops
        possible_edges = {(tf, g) for tf in tf_genes for g in all_genes if tf != g}
    else:
        possible_edges = set(permutations(unique_nodes, 2))

    all_gt_edges = set(zip(gt_no_self['Gene1'], gt_no_self['Gene2']))
    # Intersect with universe so true_edges is consistent with possible_edges
    true_edges = all_gt_edges & possible_edges

    return possible_edges, true_edges


def _compute_early_precision(
    ranked_edges: pd.DataFrame,
    true_edges: Set[Tuple[str, str]],
    possible_edges: Set[Tuple[str, str]],
    tf_edges: bool,
) -> float:
    """
    Compute the Early Precision Ratio (EPR) for one algorithm on one run.

    EPR = early_precision / random_precision, where early precision is the
    fraction of true positives among the selected top-k predictions and
    random_precision = |true_edges| / |possible_edges|. k is the number of
    true edges (the size of the ground truth network). A tie-aware cutoff
    (bestVal = max(nonZeroMin, edgeWeightTopk)) expands the selected set to
    include all tied predictions at the boundary, avoiding the large block of
    zero-weight edges when the k-th prediction has weight 0. The denominator
    is the actual number of selected predictions (not fixed k), consistent
    with the tie-expansion. Returns 0 when no predictions remain after
    filtering, and nan when the universe or ground truth is empty.

    Parameters
    ----------
    ranked_edges : pd.DataFrame
        Predicted edge list with columns Gene1, Gene2, EdgeWeight.
        Higher absolute EdgeWeight indicates greater confidence.
    true_edges : set of (str, str)
        Ground truth directed edges within the candidate universe.
    possible_edges : set of (str, str)
        Complete candidate edge universe; used to compute the random baseline
        and, when tf_edges is True, to filter predictions.
    tf_edges : bool
        When True, predictions are pre-filtered to the candidate universe
        before ranking. When False, all non-self-loop predictions are ranked.

    Returns
    -------
    float
        EPR value (>1 = better than random, =1 = at random, <1 = below random),
        0 if no predictions remain after filtering, or nan if the metric
        cannot be computed.
    """
    if not isinstance(ranked_edges, pd.DataFrame):
        raise TypeError(f"ranked_edges must be DataFrame, got {type(ranked_edges)}")
    if not isinstance(true_edges, set):
        raise TypeError(f"true_edges must be set, got {type(true_edges)}")
    if not isinstance(possible_edges, set):
        raise TypeError(f"possible_edges must be set, got {type(possible_edges)}")
    if not isinstance(tf_edges, bool):
        raise TypeError(f"tf_edges must be bool, got {type(tf_edges)}")

    num_edges = len(true_edges)
    if num_edges == 0 or len(possible_edges) == 0:
        return float('nan')

    random_eprc = num_edges / len(possible_edges)

    # Remove self-loops; deduplicate by (Gene1, Gene2) keeping highest abs weight
    predicted = ranked_edges[ranked_edges['Gene1'] != ranked_edges['Gene2']].copy()
    predicted['_abs'] = predicted['EdgeWeight'].abs()
    predicted = (
        predicted
        .sort_values('_abs', ascending=False)
        .drop_duplicates(subset=['Gene1', 'Gene2'])
        .reset_index(drop=True)
    )

    if tf_edges:
        # Restrict predictions to the TF→gene candidate universe
        keep = [(r.Gene1, r.Gene2) in possible_edges for r in predicted.itertuples()]
        predicted = predicted[keep].reset_index(drop=True)

    if predicted.empty:
        return 0.0

    # Tie-aware cutoff at rank k, capped at the number of available predictions
    maxk = min(len(predicted), num_edges)
    edge_weight_topk = float(predicted.iloc[maxk - 1]['_abs'])

    nonzero = predicted.loc[predicted['_abs'] > 0, '_abs']
    non_zero_min = float(nonzero.min()) if not nonzero.empty else 0.0

    # Expand to include all tied predictions at the boundary; avoids pulling in
    # the large block of zero-weight edges when edgeWeightTopk is 0
    best_val = max(non_zero_min, edge_weight_topk)
    selected = predicted[predicted['_abs'] >= best_val]

    if selected.empty:
        return 0.0

    selected_edges = set(zip(selected['Gene1'], selected['Gene2']))
    intersection = selected_edges & true_edges

    # Denominator is the actual selected count, not fixed k, to account for ties
    eprec = len(intersection) / len(selected_edges)
    return eprec / random_eprc


class EarlyPrecision(Evaluator):
    """
    Evaluator that computes the Early Precision Ratio (EPR) for each algorithm
    against the ground truth network.

    EPR equals early precision (fraction of true positives in top-k predictions)
    divided by the expected precision of a random predictor over the candidate
    edge universe. k equals the number of ground truth edges excluding self-loops.
    EPR > 1 indicates better-than-random performance; EPR = 1 is the random
    baseline; EPR < 1 is below random.

    Whether to restrict the candidate universe to TF→gene edges is controlled
    per dataset via the experimental_dataset field in the config (default False).
    Set experimental_dataset: true on experimental datasets and leave it unset
    (or false) for synthetic/simulated data.

    For each DatasetGroup, writes EarlyPrecision.csv to dataset_path. Rows are
    algorithms and columns are run_ids.
    """

    def __init__(self) -> None:
        pass

    def __call__(self, evaluation_data: EvaluationData) -> None:
        """
        Compute EPR per algorithm per run and write results to
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
            raise TypeError(
                f"evaluation_data must be EvaluationData, got {type(evaluation_data)}"
            )

        for dataset_group in evaluation_data:
            # Load the ground truth once per group; all runs share the same GT.
            gt_path = next(
                (run.ground_truth_path for run in dataset_group
                 if run.ground_truth_path.exists()),
                None,
            )
            if gt_path is None:
                print(
                    f"Warning: no ground truth found for dataset "
                    f"'{dataset_group.dataset_id}', skipping EarlyPrecision."
                )
                continue

            gt_df = self._load_ground_truth(gt_path)
            possible_edges, true_edges = _build_edge_universe(gt_df, dataset_group.tf_edges)

            # results[algo][run_id] = EPR score
            results: Dict[str, Dict[str, float]] = {}

            for run in dataset_group:
                for algo, ranked_edges_df in run.ranked_edges.items():
                    score = _compute_early_precision(
                        ranked_edges_df, true_edges, possible_edges, dataset_group.tf_edges
                    )
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
