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
) -> Tuple[float, float]:
    """
    Compute early precision separately for activation and inhibitory edges.

    Self-loops are excluded from predictions before ranking. For each edge
    type, the top-k predictions (by EdgeWeight) are inspected and true
    positives counted against the corresponding ground truth set. Returns
    float('nan') for a type when its k is zero or the predicted list is empty.

    Parameters
    ----------
    ranked_edges : pd.DataFrame
        Predicted edge list with columns Gene1, Gene2, EdgeWeight.
        Higher EdgeWeight indicates greater confidence.
    activation_edges : set of (str, str)
        Ground truth directed edges with Type == '+'.
    inhibitory_edges : set of (str, str)
        Ground truth directed edges with Type == '-'.
    ka : int
        Number of activation edges in the ground truth (excluding self-loops).
        Serves as both the rank cutoff and the precision denominator.
    ki : int
        Number of inhibitory edges in the ground truth (excluding self-loops).
        Serves as both the rank cutoff and the precision denominator.

    Returns
    -------
    ep_activation : float
        Fraction of true activation edges in the top-ka predictions,
        or nan if ka == 0 or no predictions remain after filtering.
    ep_inhibitory : float
        Fraction of true inhibitory edges in the top-ki predictions,
        or nan if ki == 0 or no predictions remain after filtering.
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

    # Remove self-loops — they are never meaningful predictions
    predicted = ranked_edges[ranked_edges['Gene1'] != ranked_edges['Gene2']].copy()

    if predicted.empty:
        return float('nan'), float('nan')

    # Sort once descending by edge weight; reuse for both cutoffs
    predicted_sorted = predicted.sort_values('EdgeWeight', ascending=False)

    def _ep(true_edges: Set[Tuple[str, str]], k: int) -> float:
        """Fraction of true positives in the top-k sorted predictions."""
        if k == 0:
            return float('nan')
        top_k = predicted_sorted.head(k)
        true_positives = sum(
            1 for edge in zip(top_k['Gene1'], top_k['Gene2'])
            if edge in true_edges
        )
        return true_positives / k

    ep_activation = _ep(activation_edges, ka)
    ep_inhibitory = _ep(inhibitory_edges, ki)

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

                for algo, ranked_edges_df in run.ranked_edges.items():
                    ep_act, ep_inh = _compute_signed_early_precision(
                        ranked_edges_df, activation_edges, inhibitory_edges, ka, ki
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
