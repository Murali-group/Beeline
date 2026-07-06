from itertools import product
from pathlib import Path
from typing import Dict, Iterable, Set, Tuple

import pandas as pd
from sklearn.metrics import auc, precision_recall_curve

from BLEval.data import EvaluationData
from BLEval.evaluator import Evaluator


Edge = Tuple[str, str]


def _normalize_gene_name(gene: str) -> str:
    """
    Normalize a gene name for TF-list matching.

    Matching is case-insensitive and ignores leading/trailing whitespace, while
    original gene names are preserved in candidate edges and output tables.
    """
    return str(gene).strip().upper()


def _load_tf_names(tf_list_path: Path) -> Set[str]:
    """
    Load transcription-factor names from a CSV with a required TF column.
    """
    if not isinstance(tf_list_path, Path):
        raise TypeError(f"tf_list_path must be Path, got {type(tf_list_path)}")

    if not tf_list_path.exists():
        raise FileNotFoundError(f"TF list not found: {tf_list_path}")

    tf_df = pd.read_csv(tf_list_path)

    if "TF" not in tf_df.columns:
        raise ValueError(
            f"TF list at {tf_list_path} is missing required column 'TF'"
        )

    return {
        normalized
        for normalized in (
            _normalize_gene_name(value)
            for value in tf_df["TF"].dropna()
        )
        if normalized
    }


def _prediction_lookup(
    ranked_edges: pd.DataFrame,
) -> Dict[Edge, float]:
    """
    Convert a predicted edge table into a dictionary.

    The key is the directed pair (Gene1, Gene2). The value is the absolute edge
    weight. Self-loops are removed. If an edge appears more than once, the row
    with the largest absolute weight is retained.
    """
    if not isinstance(ranked_edges, pd.DataFrame):
        raise TypeError(f"ranked_edges must be DataFrame, got {type(ranked_edges)}")

    required = {"Gene1", "Gene2", "EdgeWeight"}

    if not required.issubset(ranked_edges.columns):
        missing = sorted(required - set(ranked_edges.columns))
        raise ValueError(f"ranked_edges is missing columns: {missing}")

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

    Returns NaN when the category has only positives or only negatives, because
    the precision-recall curve is not informative in single-class cases.
    """
    if not isinstance(prediction_lookup, dict):
        raise TypeError(
            f"prediction_lookup must be dict, got {type(prediction_lookup)}"
        )
    if not isinstance(true_edges, set):
        raise TypeError(f"true_edges must be set, got {type(true_edges)}")

    candidates = sorted(set(possible_edges))

    if not candidates:
        return float("nan")

    true_labels = [
        1 if edge in true_edges else 0
        for edge in candidates
    ]

    positive_count = sum(true_labels)

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


def _construct_tf_split(
    ground_truth: pd.DataFrame,
    tf_names: Set[str],
):
    """
    Classify genes using a TF list and construct two candidate-edge universes.

    Candidate edges are restricted to edges whose source gene is classified as a
    TF. Those candidates are split into:

    1. TF source -> TF target
    2. TF source -> non-TF target

    Ground-truth self-loops are excluded from true edges. Candidate universes
    also exclude self-loops.
    """
    if not isinstance(ground_truth, pd.DataFrame):
        raise TypeError(f"ground_truth must be DataFrame, got {type(ground_truth)}")
    if not isinstance(tf_names, set):
        raise TypeError(f"tf_names must be set, got {type(tf_names)}")

    required = {"Gene1", "Gene2"}

    if not required.issubset(ground_truth.columns):
        missing = sorted(required - set(ground_truth.columns))
        raise ValueError(f"ground truth is missing columns: {missing}")

    all_edges = (
        ground_truth.loc[:, ["Gene1", "Gene2"]]
        .dropna()
        .drop_duplicates()
    )

    genes = sorted(
        set(all_edges["Gene1"]).union(set(all_edges["Gene2"]))
    )

    if len(genes) < 2:
        raise ValueError("At least two ground-truth genes are required")

    nonself_edges = all_edges.loc[
        all_edges["Gene1"] != all_edges["Gene2"]
    ]

    true_edges = set(
        zip(nonself_edges["Gene1"], nonself_edges["Gene2"])
    )

    normalized_tf_names = {
        _normalize_gene_name(name)
        for name in tf_names
    }

    tf_genes = {
        gene
        for gene in genes
        if _normalize_gene_name(gene) in normalized_tf_names
    }

    non_tf_genes = set(genes) - tf_genes

    tf_to_tf_candidates = {
        (source, target)
        for source, target in product(tf_genes, tf_genes)
        if source != target
    }

    tf_to_nontf_candidates = {
        (source, target)
        for source, target in product(tf_genes, non_tf_genes)
        if source != target
    }

    def _category_counts(candidates: Set[Edge]):
        candidate_count = len(candidates)
        positive_count = len(true_edges.intersection(candidates))
        negative_count = candidate_count - positive_count
        random_baseline = (
            positive_count / candidate_count
            if candidate_count
            else float("nan")
        )
        return (
            candidate_count,
            positive_count,
            negative_count,
            random_baseline,
        )

    (
        tf_to_tf_candidate_count,
        tf_to_tf_positive_count,
        tf_to_tf_negative_count,
        tf_to_tf_random_baseline,
    ) = _category_counts(tf_to_tf_candidates)

    (
        tf_to_nontf_candidate_count,
        tf_to_nontf_positive_count,
        tf_to_nontf_negative_count,
        tf_to_nontf_random_baseline,
    ) = _category_counts(tf_to_nontf_candidates)

    summary = {
        "GeneCount": len(genes),
        "TrueEdgeCount": len(true_edges),
        "TFNameCountInList": len(normalized_tf_names),
        "MatchedTFGeneCount": len(tf_genes),
        "MatchedTFGenes": ";".join(sorted(tf_genes)),
        "GenesClassifiedAsNonTFCount": len(non_tf_genes),
        "GenesClassifiedAsNonTF": ";".join(sorted(non_tf_genes)),
        "TFToTFCandidateCount": tf_to_tf_candidate_count,
        "TFToTFPositiveCount": tf_to_tf_positive_count,
        "TFToTFNegativeCount": tf_to_tf_negative_count,
        "TFToTFRandomBaseline": tf_to_tf_random_baseline,
        "TFToNonTFCandidateCount": tf_to_nontf_candidate_count,
        "TFToNonTFPositiveCount": tf_to_nontf_positive_count,
        "TFToNonTFNegativeCount": tf_to_nontf_negative_count,
        "TFToNonTFRandomBaseline": tf_to_nontf_random_baseline,
    }

    gene_classification = pd.DataFrame(
        {
            "Gene": genes,
            "NormalizedGene": [
                _normalize_gene_name(gene)
                for gene in genes
            ],
            "IsTF": [
                gene in tf_genes
                for gene in genes
            ],
            "ClassificationReason": [
                "matched_tf_list" if gene in tf_genes else "not_matched_tf_list"
                for gene in genes
            ],
        }
    )

    return (
        true_edges,
        tf_to_tf_candidates,
        tf_to_nontf_candidates,
        summary,
        gene_classification,
    )


class TFAUPRC(Evaluator):
    """
    Compute AUPRC separately for TF-source edges targeting TF genes and
    TF-source edges targeting non-TF genes.
    """

    def __call__(
        self,
        evaluation_data: EvaluationData,
    ) -> None:
        if not isinstance(evaluation_data, EvaluationData):
            raise TypeError(
                "evaluation_data must be an EvaluationData object"
            )

        for dataset_group in evaluation_data:
            if dataset_group.tf_list_path is None:
                print(
                    "Warning: no tf_list path configured for dataset "
                    f"'{dataset_group.dataset_id}', skipping TF AUPRC."
                )
                continue

            tf_names = _load_tf_names(dataset_group.tf_list_path)

            tf_to_tf_results: Dict[
                str,
                Dict[str, float],
            ] = {}

            tf_to_nontf_results: Dict[
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

                    split_data = _construct_tf_split(
                        ground_truth,
                        tf_names,
                    )

                (
                    true_edges,
                    tf_to_tf_candidates,
                    tf_to_nontf_candidates,
                    summary,
                    gene_classification,
                ) = split_data

                for algorithm, ranked_edges in run.ranked_edges.items():
                    lookup = _prediction_lookup(ranked_edges)

                    tf_to_tf_score = _compute_subset_auprc(
                        lookup,
                        true_edges,
                        tf_to_tf_candidates,
                    )

                    tf_to_nontf_score = _compute_subset_auprc(
                        lookup,
                        true_edges,
                        tf_to_nontf_candidates,
                    )

                    tf_to_tf_results.setdefault(
                        algorithm,
                        {},
                    )[run.run_id] = tf_to_tf_score

                    tf_to_nontf_results.setdefault(
                        algorithm,
                        {},
                    )[run.run_id] = tf_to_nontf_score

            if split_data is None:
                continue

            dataset_group.dataset_path.mkdir(
                parents=True,
                exist_ok=True,
            )

            if tf_to_tf_results:
                tf_to_tf_df = pd.DataFrame(tf_to_tf_results).T
                tf_to_tf_df.index.name = "Algorithm"
                tf_to_tf_path = (
                    dataset_group.dataset_path
                    / "TFAUPRC_TFToTF.csv"
                )
                tf_to_tf_df.to_csv(tf_to_tf_path)
                print(
                    "Wrote TF-to-TF AUPRC results to "
                    f"{tf_to_tf_path}"
                )

            if tf_to_nontf_results:
                tf_to_nontf_df = pd.DataFrame(tf_to_nontf_results).T
                tf_to_nontf_df.index.name = "Algorithm"
                tf_to_nontf_path = (
                    dataset_group.dataset_path
                    / "TFAUPRC_TFToNonTF.csv"
                )
                tf_to_nontf_df.to_csv(tf_to_nontf_path)
                print(
                    "Wrote TF-to-non-TF AUPRC results to "
                    f"{tf_to_nontf_path}"
                )

            summary_out = {
                "Dataset": dataset_group.dataset_id,
                "TFListPath": str(dataset_group.tf_list_path),
                **summary,
            }

            summary_path = (
                dataset_group.dataset_path
                / "TFAUPRC_Summary.csv"
            )
            pd.DataFrame([summary_out]).to_csv(
                summary_path,
                index=False,
            )
            print(f"Wrote TF AUPRC summary to {summary_path}")

            classification_path = (
                dataset_group.dataset_path
                / "TFAUPRC_GeneClassification.csv"
            )
            gene_classification.to_csv(
                classification_path,
                index=False,
            )
            print(
                "Wrote TF gene classification table to "
                f"{classification_path}"
            )
