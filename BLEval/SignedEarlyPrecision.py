from pathlib import Path
from typing import Dict, Set, Tuple

import pandas as pd

from BLEval.evaluator import Evaluator
from BLEval.data import EvaluationData


def _load_signed_ground_truth(
    gt_path: Path,
) -> Tuple[Set[Tuple[str, str]], Set[Tuple[str, str]]]:
    """
    Load a ground truth network and split edges by sign.

    Reads Gene1, Gene2, Type columns. Rows whose Type value is not '+' or '-'
    are silently ignored. Self-loops are retained in the returned sets; callers
    are responsible for filtering them when computing k.

    Parameters
    ----------
    gt_path : Path
        Path to a CSV with columns Gene1, Gene2, Type.
        Type must be '+' (activation) or '-' (inhibitory).

    Returns
    -------
    activation_edges : set of (str, str)
        Directed edges with Type == '+'.
    inhibitory_edges : set of (str, str)
        Directed edges with Type == '-'.
    """
    if not isinstance(gt_path, Path):
        raise TypeError(f"gt_path must be Path, got {type(gt_path)}")

    gt_df = pd.read_csv(gt_path, header=0)

    activation_edges = set(
        zip(
            gt_df.loc[gt_df['Type'] == '+', 'Gene1'],
            gt_df.loc[gt_df['Type'] == '+', 'Gene2'],
        )
    )
    inhibitory_edges = set(
        zip(
            gt_df.loc[gt_df['Type'] == '-', 'Gene1'],
            gt_df.loc[gt_df['Type'] == '-', 'Gene2'],
        )
    )

    return activation_edges, inhibitory_edges


def _compute_signed_early_precision(
    ranked_edges: pd.DataFrame,
    activation_edges: Set[Tuple[str, str]],
    inhibitory_edges: Set[Tuple[str, str]],
    ka: int,
    ki: int,
    n_possible: int,
) -> Tuple[float, float]:
    """
    Compute Early Precision Ratio separately for activation and inhibitory edges.

    Self-loops are excluded from predictions before ranking. Predictions are
    ranked by absolute EdgeWeight (magnitude = confidence). For each edge type,
    opposite-sign ground truth edges are removed from the candidate pool before
    selecting the top-k. A tie-handling cutoff (bestVal = max(nonZeroMin,
    edgeWeightTopk)) expands the selected set to include all tied predictions at
    the boundary, so the denominator is the actual number of selected predictions
    rather than a fixed k. The raw early precision is divided by the random
    predictor baseline (k / n_possible) to produce an EPR, consistent with
    EarlyPrecision.py. Returns float('nan') for a type when its k is zero,
    n_possible is zero, or no candidates remain after filtering.

    Parameters
    ----------
    ranked_edges : pd.DataFrame
        Predicted edge list with columns Gene1, Gene2, EdgeWeight.
        Higher absolute EdgeWeight indicates greater confidence.
    activation_edges : set of (str, str)
        Ground truth directed edges with Type == '+'.
    inhibitory_edges : set of (str, str)
        Ground truth directed edges with Type == '-'.
    ka : int
        Number of activation edges in the ground truth (excluding self-loops).
        Used as the rank cutoff for selecting top-k activation candidates.
    ki : int
        Number of inhibitory edges in the ground truth (excluding self-loops).
        Used as the rank cutoff for selecting top-k inhibitory candidates.
    n_possible : int
        Total number of possible directed non-self-loop edges among all genes
        in the ground truth network (n*(n-1)). Used to compute the random
        predictor baseline for each edge type.

    Returns
    -------
    ep_activation : float
        EPR for activation edges (raw precision / random baseline), or nan if
        ka == 0, n_possible == 0, or no candidates remain.
    ep_inhibitory : float
        EPR for inhibitory edges (raw precision / random baseline), or nan if
        ki == 0, n_possible == 0, or no candidates remain.
    """
    if not isinstance(ranked_edges, pd.DataFrame):
        raise TypeError(f"ranked_edges must be DataFrame, got {type(ranked_edges)}")
    if not isinstance(activation_edges, set):
        raise TypeError(f"activation_edges must be set, got {type(activation_edges)}")
    if not isinstance(inhibitory_edges, set):
        raise TypeError(f"inhibitory_edges must be set, got {type(inhibitory_edges)}")
    if not isinstance(ka, int):
        raise TypeError(f"ka must be int, got {type(ka)}")
    if not isinstance(ki, int):
        raise TypeError(f"ki must be int, got {type(ki)}")
    if not isinstance(n_possible, int):
        raise TypeError(f"n_possible must be int, got {type(n_possible)}")

    # Remove self-loops — they are never meaningful predictions
    predicted = ranked_edges[ranked_edges['Gene1'] != ranked_edges['Gene2']].copy()

    if predicted.empty:
        return float('nan'), float('nan')

    # Rank by absolute edge weight so that both positive (activating) and
    # negative (inhibitory) predictions are ordered by confidence magnitude.
    predicted['_abs'] = predicted['EdgeWeight'].abs()
    predicted_sorted = predicted.sort_values('_abs', ascending=False)

    def _ep(
        true_edges: Set[Tuple[str, str]],
        opposite_edges: Set[Tuple[str, str]],
        k: int,
    ) -> float:
        """
        Early precision for one edge type.

        Removes edges confirmed to belong to the opposite sign from the
        candidate pool, then applies a tie-aware cutoff at rank k. The
        denominator is the number of selected predictions, not a fixed k,
        so tied predictions at the boundary are counted consistently.
        """
        if k == 0:
            return float('nan')

        # Exclude predictions confirmed to be the opposite sign
        edge_tuples = list(zip(predicted_sorted['Gene1'], predicted_sorted['Gene2']))
        keep_mask = pd.Series(
            [e not in opposite_edges for e in edge_tuples],
            index=predicted_sorted.index,
        )
        candidates = predicted_sorted[keep_mask]

        if candidates.empty:
            return float('nan')

        # Cutoff at rank k; expand to include all ties at the boundary.
        # bestVal = max(nonZeroMin, edgeWeightTopk) avoids including the large
        # block of zero-weight unpredicted edges when edgeWeightTopk is 0.
        top_k = candidates.head(k)
        if top_k.empty:
            return float('nan')

        edge_weight_topk = float(top_k.iloc[-1]['_abs'])

        nonzero_weights = candidates.loc[candidates['_abs'] > 0, '_abs']
        non_zero_min = float(nonzero_weights.min()) if not nonzero_weights.empty else 0.0

        best_val = max(non_zero_min, edge_weight_topk)
        selected = candidates[candidates['_abs'] >= best_val]

        if selected.empty:
            return float('nan')

        true_positives = sum(
            1 for edge in zip(selected['Gene1'], selected['Gene2'])
            if edge in true_edges
        )
        early_precision = true_positives / len(selected)
        # Divide by random baseline to produce EPR, consistent with EarlyPrecision.py.
        random_baseline = k / n_possible
        return early_precision / random_baseline

    if n_possible == 0:
        return float('nan'), float('nan')

    ep_activation = _ep(activation_edges, inhibitory_edges, ka)
    ep_inhibitory = _ep(inhibitory_edges, activation_edges, ki)

    return ep_activation, ep_inhibitory


class SignedEarlyPrecision(Evaluator):
    """
    Evaluator that computes early precision separately for activation and
    inhibitory edges in the ground truth network.

    For activation edges, the top-ka predictions are inspected, where ka is
    the number of activation edges in the reference network excluding
    self-loops. Inhibitory precision is computed analogously with ki. For each
    DatasetGroup, writes EarlyPrecisionActivation.csv and
    EarlyPrecisionInhibitory.csv to dataset_path. Rows are algorithms and
    columns are run_ids. Runs whose ground truth file is missing or lacks a
    Type column are skipped with a warning.
    """

    def __call__(self, evaluation_data: EvaluationData) -> None:
        """
        Compute signed early precision per algorithm per run and write results
        to dataset_path/EarlyPrecisionActivation.csv and
        dataset_path/EarlyPrecisionInhibitory.csv.

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
            # results_act[algo][run_id] = activation early precision score
            # results_inh[algo][run_id] = inhibitory early precision score
            results_act: Dict[str, Dict[str, float]] = {}
            results_inh: Dict[str, Dict[str, float]] = {}

            for run in dataset_group:
                if not run.ground_truth_path.exists():
                    print(
                        f"Warning: ground truth not found at {run.ground_truth_path}, "
                        f"skipping run '{run.run_id}'."
                    )
                    continue

                # Verify Type column is present before loading signed sets
                gt_df = pd.read_csv(run.ground_truth_path, header=0, nrows=0)
                if 'Type' not in gt_df.columns:
                    print(
                        f"Warning: ground truth at {run.ground_truth_path} has no 'Type' "
                        f"column, skipping run '{run.run_id}'."
                    )
                    continue

                activation_edges, inhibitory_edges = _load_signed_ground_truth(
                    run.ground_truth_path
                )

                # k = count of each type excluding self-loops
                ka = sum(1 for g1, g2 in activation_edges if g1 != g2)
                ki = sum(1 for g1, g2 in inhibitory_edges if g1 != g2)

                # n_possible : int — total directed non-self-loop pairs among all
                # genes in the ground truth; the denominator of the random baseline.
                all_genes = set(g for e in activation_edges | inhibitory_edges for g in e)
                n = len(all_genes)
                n_possible = n * (n - 1)

                for algo, ranked_edges_df in run.ranked_edges.items():
                    ep_act, ep_inh = _compute_signed_early_precision(
                        ranked_edges_df, activation_edges, inhibitory_edges, ka, ki, n_possible
                    )
                    results_act.setdefault(algo, {})[run.run_id] = ep_act
                    results_inh.setdefault(algo, {})[run.run_id] = ep_inh

            if not results_act and not results_inh:
                continue

            dataset_group.dataset_path.mkdir(parents=True, exist_ok=True)

            # Build output DataFrames: rows = algorithms, columns = run_ids
            if results_act:
                act_df = pd.DataFrame(results_act).T
                act_df.index.name = 'Algorithm'
                act_path = dataset_group.dataset_path / 'EarlyPrecisionActivation.csv'
                act_df.to_csv(act_path)
                print(f"Wrote EarlyPrecisionActivation results to {act_path}")

            if results_inh:
                inh_df = pd.DataFrame(results_inh).T
                inh_df.index.name = 'Algorithm'
                inh_path = dataset_group.dataset_path / 'EarlyPrecisionInhibitory.csv'
                inh_df.to_csv(inh_path)
                print(f"Wrote EarlyPrecisionInhibitory results to {inh_path}")
