from typing import Dict, Set, Tuple

import networkx as nx
import pandas as pd

from BLEval.evaluator import Evaluator
from BLEval.data import EvaluationData, DatasetGroup


def _build_digraph(edges: Set[Tuple[str, str]]) -> nx.DiGraph:
    """
    Build a directed graph from a set of edges, excluding self-loops.

    Parameters
    ----------
    edges : set of (str, str)
        Directed edges as (Gene1, Gene2) tuples.

    Returns
    -------
    nx.DiGraph
    """
    if not isinstance(edges, set):
        raise TypeError(f"edges must be set, got {type(edges)}")

    G = nx.DiGraph()
    G.add_edges_from((a, b) for a, b in edges if a != b)
    return G


def _count_mi(G: nx.DiGraph) -> int:
    """
    Count two-node mutual interactions (MI) in a directed graph.

    A mutual interaction exists when both A→B and B→A are present.
    Each pair {A, B} is counted once.

    Parameters
    ----------
    G : nx.DiGraph

    Returns
    -------
    int
        Number of mutual interaction pairs.
    """
    return sum(1 for a, b in G.edges() if G.has_edge(b, a)) // 2


def _count_fbl(G: nx.DiGraph) -> int:
    """
    Count three-node feedback loops (FBL) in a directed graph.

    A feedback loop is a directed cycle of exactly three distinct nodes:
    A→B→C→A. Each cycle is counted once regardless of its starting node.

    Parameters
    ----------
    G : nx.DiGraph

    Returns
    -------
    int
        Number of three-node feedback loops.
    """
    return sum(1 for cycle in nx.simple_cycles(G) if len(cycle) == 3)


def _count_ffl(G: nx.DiGraph) -> int:
    """
    Count three-node feedforward loops (FFL) in a directed graph.

    A feedforward loop consists of three distinct nodes A, B, C where
    A→B, A→C, and B→C all exist. Node A is the source that regulates
    both B (directly and indirectly via B→C) and C directly.

    Parameters
    ----------
    G : nx.DiGraph

    Returns
    -------
    int
        Number of three-node feedforward loops.
    """
    count = 0
    for a in G.nodes():
        out_a = set(G.successors(a))
        for b in out_a:
            # C must be a successor of B, a successor of A, and distinct from A and B
            for c in G.successors(b):
                if c in out_a and c != a and c != b:
                    count += 1
    return count


def _motif_counts(edges: Set[Tuple[str, str]]) -> Tuple[int, int, int]:
    """
    Compute FBL, FFL, and MI motif counts for a set of directed edges.

    Parameters
    ----------
    edges : set of (str, str)
        Directed edges as (Gene1, Gene2) tuples.

    Returns
    -------
    tuple of (int, int, int)
        Counts of (FBL, FFL, MI) motifs.
    """
    G = _build_digraph(edges)
    return _count_fbl(G), _count_ffl(G), _count_mi(G)


def _safe_ratio(numerator: int, denominator: int) -> float:
    """
    Compute numerator / denominator, returning nan when denominator is zero.

    Parameters
    ----------
    numerator : int
    denominator : int

    Returns
    -------
    float
    """
    if denominator == 0:
        return float('nan')
    return numerator / denominator


class Motifs(Evaluator):
    """
    Evaluator that computes ratios of three-node feedback loop (FBL),
    three-node feedforward loop (FFL), and two-node mutual interaction (MI)
    motif counts between the predicted top-k network and the reference network.

    For each algorithm and run, the top-k predicted edges (k = number of
    ground truth edges, excluding self-loops) are used to build a directed
    graph. Motif counts in that graph are divided by the corresponding counts
    in the ground truth network to produce dimensionless ratios.

    For each DatasetGroup, writes three CSV files to dataset_path:
        - motifs_FBL.csv  — three-node feedback loop ratios
        - motifs_FFL.csv  — three-node feedforward loop ratios
        - motifs_MI.csv   — two-node mutual interaction ratios

    Each CSV has algorithms as rows and run_ids as columns.
    """

    def __call__(self, evaluation_data: EvaluationData) -> None:
        """
        Compute motif ratios per algorithm per run and write results to
        dataset_path/motifs_{FBL,FFL,MI}.csv.

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
            # Load the reference network from the first run with an accessible
            # ground truth. All runs in a DatasetGroup share the same ground truth.
            ref_edges: Set[Tuple[str, str]] = None
            for run in dataset_group:
                if run.ground_truth_path.exists():
                    ref_edges = self._load_ground_truth(run.ground_truth_path)
                    break

            if ref_edges is None:
                print(
                    f"Warning: no ground truth found for dataset "
                    f"'{dataset_group.dataset_id}', skipping motif computation."
                )
                continue

            # k = number of ground truth edges excluding self-loops
            k = sum(1 for g1, g2 in ref_edges if g1 != g2)
            if k == 0:
                print(
                    f"Warning: ground truth for dataset '{dataset_group.dataset_id}' "
                    "has no non-self-loop edges, skipping motif computation."
                )
                continue

            # Motif counts in the reference network
            ref_fbl, ref_ffl, ref_mi = _motif_counts(ref_edges)

            # results_X[algo][run_id] = ratio for motif X
            results_fbl: Dict[str, Dict[str, float]] = {}
            results_ffl: Dict[str, Dict[str, float]] = {}
            results_mi:  Dict[str, Dict[str, float]] = {}

            for run in dataset_group:
                for algo, ranked_edges_df in run.ranked_edges.items():
                    # Select top-k predicted edges, excluding self-loops
                    no_self = ranked_edges_df[
                        ranked_edges_df['Gene1'] != ranked_edges_df['Gene2']
                    ]
                    top_k_df = no_self.iloc[no_self['EdgeWeight'].abs().argsort()[::-1]].head(k)
                    top_k_edges = set(zip(top_k_df['Gene1'], top_k_df['Gene2']))

                    pred_fbl, pred_ffl, pred_mi = _motif_counts(top_k_edges)

                    results_fbl.setdefault(algo, {})[run.run_id] = _safe_ratio(pred_fbl, ref_fbl)
                    results_ffl.setdefault(algo, {})[run.run_id] = _safe_ratio(pred_ffl, ref_ffl)
                    results_mi.setdefault(algo,  {})[run.run_id] = _safe_ratio(pred_mi,  ref_mi)

            if not results_fbl:
                continue

            dataset_group.dataset_path.mkdir(parents=True, exist_ok=True)

            for label, results in [('FBL', results_fbl), ('FFL', results_ffl), ('MI', results_mi)]:
                out_df = pd.DataFrame(results).T
                out_df.index.name = 'Algorithm'
                out_path = dataset_group.dataset_path / f'motifs_{label}.csv'
                out_df.to_csv(out_path)
                print(f"Wrote motif ({label}) results to {out_path}")
