import math
from itertools import product
from typing import Dict, Iterable, Set, Tuple

import pandas as pd
from sklearn.metrics import auc, precision_recall_curve

from BLEval.data import EvaluationData
from BLEval.evaluator import Evaluator


Edge = Tuple[str, str]


def _prediction_lookup(
    ranked_edges: pd.DataFrame,
) -> Dict[Edge, float]:
    """
    Convert a predicted edge table into a dictionary.

    The key is the directed pair (Gene1, Gene2).
    The value is the absolute edge weight.

    Self-loops are removed. If an edge appears more than once,
    the row with the largest absolute weight is retained.
    """
    required = {"Gene1", "Gene2", "EdgeWeight"}

    if not required.issubset(ranked_edges.columns):
        missing = sorted(required - set(ranked_edges.columns))
        raise ValueError(
            f"ranked_edges is missing columns: {missing}"
        )

    predicted = ranked_edges.loc[
        ranked_edges["Gene1"] != ranked_edges["Gene2"]
    ].copy()

    predicted["_abs"] = predicted["EdgeWeight"].abs()

    predicted = (
        predicted.sort_values("_abs", ascending=False)
        .drop_duplicates(
            subset=["Gene1", "Gene2"],
            keep="first",
        )
    )

    return dict(
        zip(
            zip(predicted["Gene1"], predicted["Gene2"]),
            predicted["_abs"].astype(float),
        )
    )


def _compute_subset_auprc(
    prediction_lookup: Dict[Edge, float],
    true_edges: Set[Edge],
    possible_edges: Iterable[Edge],
) -> float:
    """
    Compute AUPRC inside one explicitly supplied candidate universe.
    """
    candidates = sorted(set(possible_edges))

    if not candidates:
        return float("nan")

    true_labels = [
        1 if edge in true_edges else 0
        for edge in candidates
    ]

    positive_count = sum(true_labels)

    # AUPRC is not informative when every candidate has the
    # same label: either no positives or no negatives.
    if positive_count == 0 or positive_count == len(candidates):
        return float("nan")

    prediction_scores = [
        prediction_lookup.get(edge, 0.0)
        for edge in candidates
    ]

    precision, recall, _ = precision_recall_curve(
        true_labels,
        prediction_scores,
    )

    return float(auc(recall, precision))


def _construct_hub_split(
    ground_truth: pd.DataFrame,
    hub_fraction: float,
):
    """
    Classify genes by ground-truth out-degree and construct two
    candidate-edge universes.

    A source gene is a hub when its out-degree is at or above the
    cutoff needed to include the requested top fraction. Ties at the
    cutoff are included.
    """
    required = {"Gene1", "Gene2"}

    if not required.issubset(ground_truth.columns):
        missing = sorted(required - set(ground_truth.columns))
        raise ValueError(
            f"ground truth is missing columns: {missing}"
        )

    if not 0 < hub_fraction <= 1:
        raise ValueError(
            "hub_fraction must be greater than 0 and at most 1"
        )

    edges = (
        ground_truth.loc[
            ground_truth["Gene1"] != ground_truth["Gene2"],
            ["Gene1", "Gene2"],
        ]
        .dropna()
        .drop_duplicates()
    )

    genes = sorted(
        set(edges["Gene1"]).union(set(edges["Gene2"]))
    )

    if len(genes) < 2:
        raise ValueError(
            "At least two ground-truth genes are required"
        )

    true_edges = set(zip(edges["Gene1"], edges["Gene2"]))

    out_degree = (
        edges.groupby("Gene1")
        .size()
        .reindex(genes, fill_value=0)
        .astype(int)
    )

    ranked_degrees = out_degree.sort_values(
        ascending=False,
        kind="stable",
    )

    requested_hub_count = max(
        1,
        math.ceil(hub_fraction * len(genes)),
    )

    cutoff = int(
        ranked_degrees.iloc[requested_hub_count - 1]
    )

    hub_genes = set(
        out_degree[out_degree >= cutoff].index
    )

    nonhub_genes = set(genes) - hub_genes

    hub_candidates = {
        (source, target)
        for source, target in product(hub_genes, genes)
        if source != target
    }

    nonhub_candidates = {
        (source, target)
        for source, target in product(nonhub_genes, genes)
        if source != target
    }

    hub_positive_count = len(
        true_edges.intersection(hub_candidates)
    )

    nonhub_positive_count = len(
        true_edges.intersection(nonhub_candidates)
    )

    hub_candidate_count = len(hub_candidates)
    nonhub_candidate_count = len(nonhub_candidates)

    summary = {
        "HubFraction": hub_fraction,
        "GeneCount": len(genes),
        "TrueEdgeCount": len(true_edges),
        "RequestedHubCount": requested_hub_count,
        "OutDegreeCutoff": cutoff,
        "HubCountAfterTies": len(hub_genes),
        "HubGenes": ";".join(sorted(hub_genes)),
        "NonHubCount": len(nonhub_genes),
        "HubCandidateCount": hub_candidate_count,
        "HubPositiveCount": hub_positive_count,
        "HubNegativeCount": (
            hub_candidate_count - hub_positive_count
        ),
        "HubRandomBaseline": (
            hub_positive_count / hub_candidate_count
            if hub_candidate_count
            else float("nan")
        ),
        "NonHubCandidateCount": nonhub_candidate_count,
        "NonHubPositiveCount": nonhub_positive_count,
        "NonHubNegativeCount": (
            nonhub_candidate_count
            - nonhub_positive_count
        ),
        "NonHubRandomBaseline": (
            nonhub_positive_count
            / nonhub_candidate_count
            if nonhub_candidate_count
            else float("nan")
        ),
    }

    degree_table = (
        out_degree
        .sort_values(ascending=False, kind="stable")
        .rename_axis("Gene")
        .reset_index(name="OutDegree")
    )

    degree_table["IsHub"] = degree_table["Gene"].isin(
        hub_genes
    )

    return (
        true_edges,
        hub_candidates,
        nonhub_candidates,
        summary,
        degree_table,
    )


class HubAUPRC(Evaluator):
    """
    Compute AUPRC separately for edges whose source is a hub and
    edges whose source is not a hub.
    """

    def __init__(self, hub_fraction: float = 0.20) -> None:
        if not 0 < hub_fraction <= 1:
            raise ValueError(
                "hub_fraction must be greater than 0 and at most 1"
            )

        self.hub_fraction = hub_fraction

    def __call__(
        self,
        evaluation_data: EvaluationData,
    ) -> None:
        if not isinstance(evaluation_data, EvaluationData):
            raise TypeError(
                "evaluation_data must be an EvaluationData object"
            )

        fraction_label = f"{self.hub_fraction:.2f}"

        for dataset_group in evaluation_data:
            hub_results: Dict[
                str,
                Dict[str, float],
            ] = {}

            nonhub_results: Dict[
                str,
                Dict[str, float],
            ] = {}

            split_data = None

            for run in dataset_group:
                if not run.ground_truth_path.exists():
                    print(
                        "Warning: ground truth not found at "
                        f"{run.ground_truth_path}; "
                        f"skipping run '{run.run_id}'."
                    )
                    continue

                if split_data is None:
                    ground_truth = self._load_ground_truth(
                        run.ground_truth_path
                    )

                    split_data = _construct_hub_split(
                        ground_truth,
                        self.hub_fraction,
                    )

                (
                    true_edges,
                    hub_candidates,
                    nonhub_candidates,
                    summary,
                    degree_table,
                ) = split_data

                for algorithm, ranked_edges in (
                    run.ranked_edges.items()
                ):
                    lookup = _prediction_lookup(ranked_edges)

                    hub_score = _compute_subset_auprc(
                        lookup,
                        true_edges,
                        hub_candidates,
                    )

                    nonhub_score = _compute_subset_auprc(
                        lookup,
                        true_edges,
                        nonhub_candidates,
                    )

                    hub_results.setdefault(
                        algorithm,
                        {},
                    )[run.run_id] = hub_score

                    nonhub_results.setdefault(
                        algorithm,
                        {},
                    )[run.run_id] = nonhub_score

            if split_data is None:
                continue

            dataset_group.dataset_path.mkdir(
                parents=True,
                exist_ok=True,
            )

            if hub_results:
                hub_df = pd.DataFrame(hub_results).T
                hub_df.index.name = "Algorithm"

                hub_path = (
                    dataset_group.dataset_path
                    / f"AUPRCHub-f{fraction_label}.csv"
                )

                hub_df.to_csv(hub_path)
                print(
                    f"Wrote hub AUPRC results to {hub_path}"
                )

            if nonhub_results:
                nonhub_df = pd.DataFrame(
                    nonhub_results
                ).T

                nonhub_df.index.name = "Algorithm"

                nonhub_path = (
                    dataset_group.dataset_path
                    / f"AUPRCNonHub-f{fraction_label}.csv"
                )

                nonhub_df.to_csv(nonhub_path)
                print(
                    "Wrote non-hub AUPRC results to "
                    f"{nonhub_path}"
                )

            summary_path = (
                dataset_group.dataset_path
                / f"HubSummary-f{fraction_label}.csv"
            )

            pd.DataFrame([summary]).to_csv(
                summary_path,
                index=False,
            )

            print(
                f"Wrote hub summary to {summary_path}"
            )

            degree_path = (
                dataset_group.dataset_path
                / f"HubDegrees-f{fraction_label}.csv"
            )

            degree_table.to_csv(
                degree_path,
                index=False,
            )

            print(
                f"Wrote hub degree table to {degree_path}"
            )
