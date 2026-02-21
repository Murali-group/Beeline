from typing import Dict, Set, Tuple

import pandas as pd
from sklearn.metrics import auc, roc_curve

from BLEval.evaluator import Evaluator
from BLEval.data import EvaluationData


def _compute_auroc(
    ranked_edges: pd.DataFrame,
    true_edges: Set[Tuple[str, str]],
) -> float:
    """
    Compute the area under the ROC curve for one algorithm on one run.

    Self-loops are excluded before scoring. Returns float('nan') if the
    ROC curve is undefined (no positive or no negative examples in the
    predicted list).

    Parameters
    ----------
    ranked_edges : pd.DataFrame
        Predicted edge list with columns Gene1, Gene2, EdgeWeight.
        Higher EdgeWeight indicates greater confidence.
    true_edges : set of (str, str)
        Ground truth directed edges.

    Returns
    -------
    float
        AUROC score in [0, 1], or nan if the metric cannot be computed.
    """
    if not isinstance(ranked_edges, pd.DataFrame):
        raise TypeError(f"ranked_edges must be DataFrame, got {type(ranked_edges)}")
    if not isinstance(true_edges, set):
        raise TypeError(f"true_edges must be set, got {type(true_edges)}")

    # Remove self-loops — they are never meaningful predictions
    predicted = ranked_edges[ranked_edges['Gene1'] != ranked_edges['Gene2']].copy()

    if predicted.empty:
        return float('nan')

    labels = [
        1 if edge in true_edges else 0
        for edge in zip(predicted['Gene1'], predicted['Gene2'])
    ]
    scores = predicted['EdgeWeight'].values

    # ROC curve is undefined when only one class appears in the predictions
    if sum(labels) == 0 or sum(labels) == len(labels):
        return float('nan')

    fpr, tpr, _ = roc_curve(labels, scores)
    return float(auc(fpr, tpr))


class AUROC(Evaluator):
    """
    Evaluator that computes the area under the ROC curve (AUROC)
    for each algorithm against the ground truth network.

    For each DatasetGroup, writes AUROC.csv to dataset_path. Rows are
    algorithms and columns are run_ids. Runs whose ground truth file is
    missing are skipped with a warning.
    """

    def __call__(self, evaluation_data: EvaluationData) -> None:
        """
        Compute AUROC per algorithm per run and write results to dataset_path/AUROC.csv.

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
            # results[algo][run_id] = auroc score
            results: Dict[str, Dict[str, float]] = {}

            for run in dataset_group:
                if not run.ground_truth_path.exists():
                    print(
                        f"Warning: ground truth not found at {run.ground_truth_path}, "
                        f"skipping run '{run.run_id}'."
                    )
                    continue

                true_edges = self._load_ground_truth(run.ground_truth_path)

                for algo, ranked_edges_df in run.ranked_edges.items():
                    score = _compute_auroc(ranked_edges_df, true_edges)
                    results.setdefault(algo, {})[run.run_id] = score

            if not results:
                continue

            # Build output DataFrame: rows = algorithms, columns = run_ids
            out_df = pd.DataFrame(results).T
            out_df.index.name = 'Algorithm'

            dataset_group.dataset_path.mkdir(parents=True, exist_ok=True)
            out_path = dataset_group.dataset_path / 'AUROC.csv'
            out_df.to_csv(out_path)
            print(f"Wrote AUROC results to {out_path}")
