from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from BLEval.evaluator import Evaluator
from BLEval.data import EvaluationData


def _compute_borda_aggregations(
    ranked_edges: Dict[str, pd.DataFrame],
) -> Dict[str, pd.DataFrame]:
    """
    Compute four Borda count aggregations over a set of algorithm predictions.

    For each algorithm, every predicted edge (excluding self-loops) is given a
    rank (1 = highest EdgeWeight). Edges absent from an algorithm's list are
    assigned a rank one beyond that algorithm's last predicted edge. The Borda
    score for an edge from algorithm i is:

        borda_i(edge) = N_total - rank_i(edge) + 1

    where N_total is the total number of unique non-self-loop edges across all
    algorithms. Scores are then aggregated with mean, median, min, and max
    across algorithms. An empty dict is returned when there are no predictions.

    Parameters
    ----------
    ranked_edges : dict[str, pd.DataFrame]
        Keyed by algorithm name. Each DataFrame has columns Gene1, Gene2,
        EdgeWeight. Higher EdgeWeight means greater confidence.

    Returns
    -------
    dict[str, pd.DataFrame]
        Keyed by aggregation method name: 'mean', 'median', 'min', 'max'.
        Each value is a DataFrame with columns Gene1, Gene2, Score, sorted by
        Score descending.
    """
    if not isinstance(ranked_edges, dict):
        raise TypeError(f"ranked_edges must be dict, got {type(ranked_edges)}")

    if not ranked_edges:
        return {}

    # Collect the union of all non-self-loop edges across all algorithms
    all_edges: set = set()
    for df in ranked_edges.values():
        if not isinstance(df, pd.DataFrame):
            raise TypeError(f"Each ranked_edges value must be DataFrame, got {type(df)}")
        no_loops = df[df['Gene1'] != df['Gene2']]
        all_edges.update(zip(no_loops['Gene1'], no_loops['Gene2']))

    if not all_edges:
        return {}

    # Deterministic ordering so the score matrix rows are stable
    edge_list: List[Tuple[str, str]] = sorted(all_edges)
    n_total = len(edge_list)
    edge_index = {edge: idx for idx, edge in enumerate(edge_list)}

    # scores_matrix[edge_idx, algo_idx] = Borda score for that edge from that algo
    n_algos = len(ranked_edges)
    scores_matrix = np.zeros((n_total, n_algos), dtype=float)

    for col_idx, (algo, df) in enumerate(ranked_edges.items()):
        no_loops = df[df['Gene1'] != df['Gene2']].copy()
        # Rank edges: sort descending by EdgeWeight, assign 1-based rank
        no_loops = no_loops.sort_values('EdgeWeight', ascending=False).reset_index(drop=True)
        n_algo = len(no_loops)

        # Borda score for predicted edges (rank = position index + 1)
        for rank_minus1, row in enumerate(no_loops.itertuples(index=False)):
            edge = (row.Gene1, row.Gene2)
            borda = n_total - (rank_minus1 + 1) + 1  # = n_total - rank_minus1
            scores_matrix[edge_index[edge], col_idx] = borda

        # Borda score for edges absent from this algorithm's list:
        # treat their rank as n_algo + 1 → score = n_total - n_algo
        absent_score = n_total - n_algo
        predicted_set = set(zip(no_loops['Gene1'], no_loops['Gene2']))
        for edge in edge_list:
            if edge not in predicted_set:
                scores_matrix[edge_index[edge], col_idx] = absent_score

    gene1 = [e[0] for e in edge_list]
    gene2 = [e[1] for e in edge_list]

    aggregations = {
        'mean':   np.mean(scores_matrix,   axis=1),
        'median': np.median(scores_matrix, axis=1),
        'min':    np.min(scores_matrix,    axis=1),
        'max':    np.max(scores_matrix,    axis=1),
    }

    results: Dict[str, pd.DataFrame] = {}
    for method, agg_scores in aggregations.items():
        out_df = pd.DataFrame({'Gene1': gene1, 'Gene2': gene2, 'Score': agg_scores})
        out_df = out_df.sort_values('Score', ascending=False).reset_index(drop=True)
        results[method] = out_df

    return results


class Borda(Evaluator):
    """
    Evaluator that aggregates per-algorithm ranked edge lists into a single
    consensus ranking using Borda count.

    For each run, four aggregation methods are computed: mean, median, min, and
    max of per-algorithm Borda scores. Results across all runs in a DatasetGroup
    are combined into a single dataset-level DataFrame per method, where rows
    are edges (Gene1|Gene2) and columns are run_ids. Output is written to
    dataset_path as:

        dataset_path/BordaMean.csv
        dataset_path/BordaMedian.csv
        dataset_path/BordaMin.csv
        dataset_path/BordaMax.csv

    Each file has a row for every unique edge seen across all runs, with the
    aggregated Borda score for each run as a column (NaN where an edge does not
    appear in a particular run). Runs with no algorithm predictions are skipped.
    """

    _METHODS = ('mean', 'median', 'min', 'max')

    def __call__(self, evaluation_data: EvaluationData) -> None:
        """
        Compute Borda aggregations per run, combine across runs, and write one
        CSV per aggregation method to each DatasetGroup's dataset_path.

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
            # method → {run_id → Series(Score, index=(Gene1, Gene2))}
            per_method: Dict[str, Dict[str, pd.Series]] = {m: {} for m in self._METHODS}

            for run in dataset_group:
                if not run.ranked_edges:
                    print(
                        f"Warning: no algorithm predictions for run '{run.run_id}' "
                        f"in dataset '{dataset_group.dataset_id}', skipping."
                    )
                    continue

                aggregations = _compute_borda_aggregations(run.ranked_edges)

                if not aggregations:
                    print(
                        f"Warning: Borda aggregation produced no output for run "
                        f"'{run.run_id}' in dataset '{dataset_group.dataset_id}', skipping."
                    )
                    continue

                for method, ranked_df in aggregations.items():
                    # Index the Score series by (Gene1, Gene2) tuples for alignment
                    idx = pd.MultiIndex.from_arrays(
                        [ranked_df['Gene1'], ranked_df['Gene2']], names=['Gene1', 'Gene2']
                    )
                    per_method[method][run.run_id] = pd.Series(
                        ranked_df['Score'].values, index=idx, name=run.run_id
                    )

            # Skip dataset if no runs produced any output
            if not any(per_method[m] for m in self._METHODS):
                continue

            dataset_group.dataset_path.mkdir(parents=True, exist_ok=True)

            for method in self._METHODS:
                run_series = per_method[method]
                if not run_series:
                    continue

                # Combine runs: rows = (Gene1, Gene2), columns = run_ids
                # Missing edges for a run are filled with NaN
                out_df = pd.concat(run_series.values(), axis=1)
                out_df.index.names = ['Gene1', 'Gene2']

                filename = f'Borda{method.capitalize()}.csv'
                out_path = dataset_group.dataset_path / filename
                out_df.to_csv(out_path)
                print(f"Wrote {filename} to {out_path}")
