from itertools import permutations
from typing import Dict

import pandas as pd
from sklearn.metrics import auc, roc_curve

from BLEval.evaluator import Evaluator
from BLEval.data import EvaluationData


def _compute_auroc(
    ranked_edges: pd.DataFrame,
    gt_df: pd.DataFrame,
) -> float:
    """
    Compute the area under the ROC curve for one algorithm on one run.

    Scores are computed over the fixed universe of all possible directed pairs
    among ground truth genes (all permutations of length 2, self-loops excluded).
    Predicted edges absent from the universe receive score 0. Edge weights are
    taken as absolute values so that both positive (activating) and negative
    (inhibitory) predictions are ranked by confidence magnitude. Duplicate
    (Gene1, Gene2) pairs in ranked_edges are deduplicated by keeping the row with
    the highest absolute EdgeWeight. Returns float('nan') if the ROC curve is
    undefined (only one class present in the universe).

    Parameters
    ----------
    ranked_edges : pd.DataFrame
        Predicted edge list with columns Gene1, Gene2, EdgeWeight.
        Higher absolute EdgeWeight indicates greater confidence.
    gt_df : pd.DataFrame
        Ground truth network with columns Gene1, Gene2. Used to derive both the
        set of true edges and the fixed scoring universe.

    Returns
    -------
    float
        AUROC score in [0, 1], or nan if the metric cannot be computed.
    """
    if not isinstance(ranked_edges, pd.DataFrame):
        raise TypeError(f"ranked_edges must be DataFrame, got {type(ranked_edges)}")
    if not isinstance(gt_df, pd.DataFrame):
        raise TypeError(f"gt_df must be DataFrame, got {type(gt_df)}")

    # Fixed universe: all directed pairs among GT genes, self-loops excluded
    gt_genes = sorted(set(gt_df['Gene1']).union(set(gt_df['Gene2'])))
    possible_edges = list(permutations(gt_genes, 2))

    true_edge_set = set(zip(gt_df['Gene1'], gt_df['Gene2']))

    # Deduplicate predicted edges: per (Gene1, Gene2), keep highest abs weight.
    # Self-loops removed before deduplication as they are never valid predictions.
    pred = ranked_edges[ranked_edges['Gene1'] != ranked_edges['Gene2']].copy()
    pred['_abs'] = pred['EdgeWeight'].abs()
    pred = pred.sort_values('_abs', ascending=False).drop_duplicates(
        subset=['Gene1', 'Gene2']
    )
    pred_lookup = dict(zip(zip(pred['Gene1'], pred['Gene2']), pred['_abs'].astype(float)))

    true_labels = [1 if e in true_edge_set else 0 for e in possible_edges]
    pred_scores = [pred_lookup.get(e, 0.0) for e in possible_edges]

    # ROC curve is undefined when only one class appears in the universe
    n_pos = sum(true_labels)
    if n_pos == 0 or n_pos == len(true_labels):
        return float('nan')

    fpr, tpr, _ = roc_curve(true_labels, pred_scores)
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

                gt_df = self._load_ground_truth(run.ground_truth_path)

                for algo, ranked_edges_df in run.ranked_edges.items():
                    score = _compute_auroc(ranked_edges_df, gt_df)
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
