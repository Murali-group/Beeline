import math
import tempfile
import unittest
from itertools import product
from pathlib import Path

import pandas as pd

from BLEval.TFAUPRC import (
    _compute_subset_auprc,
    _construct_tf_split,
    _load_tf_names,
    _normalize_gene_name,
    _prediction_lookup,
)


class TestTFAUPRC(unittest.TestCase):
    def test_normalize_gene_name_is_case_insensitive_and_strips_whitespace(self):
        self.assertEqual(_normalize_gene_name(" sox9 "), "SOX9")
        self.assertEqual(_normalize_gene_name("Sox9"), "SOX9")
        self.assertEqual(_normalize_gene_name("SOX9"), "SOX9")

    def test_load_tf_names_requires_tf_column_and_matches_case_insensitively(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tf_path = Path(tmpdir) / "tfs.csv"
            pd.DataFrame(
                {
                    "TF": [" a ", "B", "c"],
                    "Sources": ["x", "y", "z"],
                }
            ).to_csv(tf_path, index=False)

            tf_names = _load_tf_names(tf_path)

        self.assertEqual(tf_names, {"A", "B", "C"})

    def test_load_tf_names_missing_file_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            missing_path = Path(tmpdir) / "missing.csv"

            with self.assertRaises(FileNotFoundError):
                _load_tf_names(missing_path)

    def test_load_tf_names_malformed_csv_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tf_path = Path(tmpdir) / "bad_tfs.csv"
            pd.DataFrame(
                {
                    "Gene": ["A", "B"],
                    "Sources": ["x", "y"],
                }
            ).to_csv(tf_path, index=False)

            with self.assertRaises(ValueError):
                _load_tf_names(tf_path)

    def test_tf_split_partitions_tf_source_candidate_edges(self):
        ground_truth = pd.DataFrame(
            {
                "Gene1": ["A", "A", "B", "C"],
                "Gene2": ["B", "C", "C", "A"],
            }
        )

        # A and B are TFs. C is not.
        tf_names = {"a", "b"}

        (
            true_edges,
            tf_to_tf_candidates,
            tf_to_nontf_candidates,
            summary,
            gene_classification,
        ) = _construct_tf_split(ground_truth, tf_names)

        genes = {"A", "B", "C"}
        tf_genes = {"A", "B"}
        non_tf_genes = {"C"}

        expected_tf_to_tf = {
            (source, target)
            for source, target in product(tf_genes, tf_genes)
            if source != target
        }

        expected_tf_to_nontf = {
            (source, target)
            for source, target in product(tf_genes, non_tf_genes)
            if source != target
        }

        expected_all_tf_source = {
            (source, target)
            for source, target in product(tf_genes, genes)
            if source != target
        }

        self.assertEqual(tf_to_tf_candidates, expected_tf_to_tf)
        self.assertEqual(tf_to_nontf_candidates, expected_tf_to_nontf)

        self.assertTrue(
            tf_to_tf_candidates.isdisjoint(tf_to_nontf_candidates)
        )

        self.assertEqual(
            tf_to_tf_candidates | tf_to_nontf_candidates,
            expected_all_tf_source,
        )

        self.assertTrue(
            all(source in tf_genes for source, _ in tf_to_tf_candidates)
        )
        self.assertTrue(
            all(source in tf_genes for source, _ in tf_to_nontf_candidates)
        )

        self.assertTrue(
            all(target in tf_genes for _, target in tf_to_tf_candidates)
        )
        self.assertTrue(
            all(target in non_tf_genes for _, target in tf_to_nontf_candidates)
        )

        self.assertNotIn(("C", "A"), tf_to_tf_candidates)
        self.assertNotIn(("C", "A"), tf_to_nontf_candidates)

        self.assertEqual(true_edges, {("A", "B"), ("A", "C"), ("B", "C"), ("C", "A")})
        self.assertEqual(summary["MatchedTFGeneCount"], 2)
        self.assertEqual(summary["GenesClassifiedAsNonTFCount"], 1)
        self.assertEqual(set(gene_classification["Gene"]), genes)

    def test_tf_split_removes_self_loops_from_true_edges_and_candidates(self):
        ground_truth = pd.DataFrame(
            {
                "Gene1": ["A", "A", "B"],
                "Gene2": ["A", "B", "A"],
            }
        )

        (
            true_edges,
            tf_to_tf_candidates,
            tf_to_nontf_candidates,
            summary,
            gene_classification,
        ) = _construct_tf_split(ground_truth, {"A", "B"})

        self.assertNotIn(("A", "A"), true_edges)
        self.assertNotIn(("A", "A"), tf_to_tf_candidates)
        self.assertNotIn(("B", "B"), tf_to_tf_candidates)
        self.assertEqual(true_edges, {("A", "B"), ("B", "A")})

    def test_prediction_lookup_removes_self_loops_and_keeps_largest_abs_duplicate(self):
        ranked_edges = pd.DataFrame(
            {
                "Gene1": ["A", "A", "A", "B"],
                "Gene2": ["B", "B", "A", "A"],
                "EdgeWeight": [0.2, -0.9, 5.0, -0.3],
            }
        )

        lookup = _prediction_lookup(ranked_edges)

        self.assertNotIn(("A", "A"), lookup)
        self.assertEqual(lookup[("A", "B")], 0.9)
        self.assertEqual(lookup[("B", "A")], 0.3)

    def test_prediction_lookup_requires_columns(self):
        ranked_edges = pd.DataFrame(
            {
                "Gene1": ["A"],
                "Gene2": ["B"],
            }
        )

        with self.assertRaises(ValueError):
            _prediction_lookup(ranked_edges)

    def test_ground_truth_requires_gene_columns(self):
        ground_truth = pd.DataFrame(
            {
                "Source": ["A"],
                "Target": ["B"],
            }
        )

        with self.assertRaises(ValueError):
            _construct_tf_split(ground_truth, {"A"})

    def test_no_positive_edges_returns_nan(self):
        score = _compute_subset_auprc(
            prediction_lookup={
                ("A", "B"): 0.8,
                ("A", "C"): 0.4,
            },
            true_edges=set(),
            possible_edges={
                ("A", "B"),
                ("A", "C"),
            },
        )

        self.assertTrue(math.isnan(score))

    def test_no_negative_edges_returns_nan(self):
        candidates = {
            ("A", "B"),
            ("A", "C"),
        }

        score = _compute_subset_auprc(
            prediction_lookup={
                ("A", "B"): 0.8,
                ("A", "C"): 0.4,
            },
            true_edges=candidates,
            possible_edges=candidates,
        )

        self.assertTrue(math.isnan(score))

    def test_mixed_labels_produce_valid_auprc(self):
        score = _compute_subset_auprc(
            prediction_lookup={
                ("A", "B"): 0.9,
                ("A", "C"): 0.7,
                ("B", "A"): 0.2,
                ("B", "C"): 0.1,
            },
            true_edges={
                ("A", "B"),
                ("B", "C"),
            },
            possible_edges={
                ("A", "B"),
                ("A", "C"),
                ("B", "A"),
                ("B", "C"),
            },
        )

        self.assertFalse(math.isnan(score))
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)


if __name__ == "__main__":
    unittest.main()
