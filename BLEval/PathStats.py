from typing import Dict, Set, Tuple

import networkx as nx
import pandas as pd

from BLEval.evaluator import Evaluator
from BLEval.data import EvaluationData


def _build_ref_graph(true_edges: Set[Tuple[str, str]]) -> nx.DiGraph:
    """
    Build a directed reference graph from a set of ground truth edges.

    Self-loops are excluded since they are not meaningful for path analysis.

    Parameters
    ----------
    true_edges : set of (str, str)
        Ground truth directed edges.

    Returns
    -------
    nx.DiGraph
        Directed graph containing all non-self-loop ground truth edges.
    """
    if not isinstance(true_edges, set):
        raise TypeError(f"true_edges must be set, got {type(true_edges)}")

    G = nx.DiGraph()
    G.add_edges_from((u, v) for u, v in true_edges if u != v)
    return G


def _path_stats(
    pred_graph: nx.DiGraph,
    ref_graph: nx.DiGraph,
) -> Dict:
    """
    Compute path statistics between a predicted graph and a reference graph.

    For each false positive edge (u, v) in the predicted graph, the shortest
    directed path from u to v in the reference graph is looked up. Path lengths
    are binned into {2, 3, 4, 5}; edges with no path in the reference graph
    increment numFP_noPath. The key 0 is included for compatibility with the
    reference format but is never incremented (path length 0 implies u == v,
    which is excluded as a self-loop).

    Parameters
    ----------
    pred_graph : nx.DiGraph
        Top-k predicted network.
    ref_graph : nx.DiGraph
        Ground truth directed network.

    Returns
    -------
    dict
        Keys: '0', '2', '3', '4', '5' (path length counts), numPred, numTP,
        numFP_withPath, numFP_noPath.
    """
    if not isinstance(pred_graph, nx.DiGraph):
        raise TypeError(f"pred_graph must be nx.DiGraph, got {type(pred_graph)}")
    if not isinstance(ref_graph, nx.DiGraph):
        raise TypeError(f"ref_graph must be nx.DiGraph, got {type(ref_graph)}")

    ref_edges  = set(ref_graph.edges())
    pred_edges = set(pred_graph.edges())

    true_positives  = pred_edges & ref_edges
    false_positives = pred_edges - ref_edges

    nopath   = 0
    yespath  = 0
    path_length_counts = {0: 0, 2: 0, 3: 0, 4: 0, 5: 0}

    for u, v in false_positives:
        try:
            path = nx.dijkstra_path(ref_graph, u, v)
            path_len = len(path) - 1
            yespath += 1
            if path_len in path_length_counts:
                path_length_counts[path_len] += 1
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            nopath += 1

    return {
        **{str(k): v for k, v in path_length_counts.items()},
        'numPred':        len(pred_edges),
        'numTP':          len(true_positives),
        'numFP_withPath': yespath,
        'numFP_noPath':   nopath,
    }


def _compute_path_stats(
    ranked_edges: pd.DataFrame,
    ref_graph: nx.DiGraph,
) -> Dict:
    """
    Compute path statistics for one algorithm on one run.

    Filters self-loops, zero-weight edges, and duplicates before taking the
    top-k predictions (k = number of reference graph edges). Returns None when
    there are no predictions or the reference graph is empty.

    Parameters
    ----------
    ranked_edges : pd.DataFrame
        Predicted edge list with columns Gene1, Gene2, EdgeWeight.
    ref_graph : nx.DiGraph
        Ground truth directed network used as the path reference.

    Returns
    -------
    dict or None
        Path statistics dict (see _path_stats), or None if the metric cannot
        be computed.
    """
    if not isinstance(ranked_edges, pd.DataFrame):
        raise TypeError(f"ranked_edges must be DataFrame, got {type(ranked_edges)}")
    if not isinstance(ref_graph, nx.DiGraph):
        raise TypeError(f"ref_graph must be nx.DiGraph, got {type(ref_graph)}")

    k = len(ref_graph.edges())
    if k == 0:
        return None

    # Filter self-loops, round and take absolute value, drop zeros and duplicates
    pred = ranked_edges[ranked_edges['Gene1'] != ranked_edges['Gene2']].copy()
    pred['EdgeWeight'] = pred['EdgeWeight'].abs().round(6)
    pred = pred[pred['EdgeWeight'] > 0]
    pred.drop_duplicates(keep='first', inplace=True)
    pred.reset_index(drop=True, inplace=True)

    if pred.empty:
        return None

    # Take top-k by EdgeWeight; use the weight at position k-1 as the cutoff
    maxk      = min(len(pred), k)
    threshold = pred.iloc[maxk - 1]['EdgeWeight']
    top_k     = pred[pred['EdgeWeight'] >= threshold]

    pred_graph = nx.DiGraph()
    pred_graph.add_edges_from(zip(top_k['Gene1'], top_k['Gene2']))

    return _path_stats(pred_graph, ref_graph)


class PathStats(Evaluator):
    """
    Evaluator that characterises false positive predicted edges by their
    shortest-path distance in the ground truth network.

    For each algorithm, the top-k predicted edges (k = |reference edges|) are
    split into true positives and false positives. Each false positive (u, v)
    is looked up in the reference graph; if a directed path from u to v exists
    its length is binned (2–5); otherwise it is counted as having no path.

    For each DatasetGroup, one CSV per run is written to dataset_path:

        dataset_path/PathStats_{run_id}.csv

    Rows are algorithms; columns are '0', '2', '3', '4', '5', numPred, numTP,
    numFP_withPath, numFP_noPath. Runs whose ground truth file is missing are
    skipped with a warning.
    """

    def __call__(self, evaluation_data: EvaluationData) -> None:
        """
        Compute path statistics per algorithm per run and write results to
        dataset_path/PathStats_{run_id}.csv.

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
            dataset_group.dataset_path.mkdir(parents=True, exist_ok=True)

            for run in dataset_group:
                if not run.ground_truth_path.exists():
                    print(
                        f"Warning: ground truth not found at {run.ground_truth_path}, "
                        f"skipping run '{run.run_id}'."
                    )
                    continue

                true_edges = self._load_ground_truth(run.ground_truth_path)
                ref_graph  = _build_ref_graph(true_edges)

                # results[algo] = path stats dict
                results: Dict[str, Dict] = {}
                for algo, ranked_edges_df in run.ranked_edges.items():
                    stats = _compute_path_stats(ranked_edges_df, ref_graph)
                    if stats is not None:
                        results[algo] = stats

                if not results:
                    continue

                out_df = pd.DataFrame(results).T
                out_df.index.name = 'Algorithm'

                out_path = dataset_group.dataset_path / f'PathStats_{run.run_id}.csv'
                out_df.to_csv(out_path)
                print(f"Wrote PathStats results to {out_path}")
