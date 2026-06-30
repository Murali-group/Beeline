import math
import unittest
from itertools import product

import pandas as pd

from BLEval.HubAUPRC import (
    _compute_subset_auprc,
    _construct_hub_split,
    _prediction_lookup,
)


class TestHubAUPRC(unittest.TestCase):
    def test_prediction_lookup_removes_self_loops_and_duplicates(self):
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

    def test_hub_split_includes_ties_and_partitions_all_pairs(self):
        ground_truth = pd.DataFrame(
            {
                "Gene1": ["A", "A", "A", "B", "C"],
                "Gene2": ["B", "C", "D", "C", "D"],
            }
        )

        (
            true_edges,
            hub_candidates,
            nonhub_candidates,
            summary,
            degree_table,
        ) = _construct_hub_split(
            ground_truth,
            hub_fraction=0.50,
        )

        # Out-degrees are A=3, B=1, C=1, D=0.
        # Two hubs are requested, but B and C tie at the cutoff,
        # so three genes are included.
        self.assertEqual(summary["RequestedHubCount"], 2)
        self.assertEqual(summary["OutDegreeCutoff"], 1)
        self.assertEqual(summary["HubCountAfterTies"], 3)
        self.assertEqual(summary["HubGenes"], "A;B;C")

        genes = {"A", "B", "C", "D"}

        all_nonself_pairs = {
            (source, target)
            for source, target in product(genes, genes)
            if source != target
        }

        self.assertTrue(
            hub_candidates.isdisjoint(nonhub_candidates)
        )

        self.assertEqual(
            hub_candidates | nonhub_candidates,
            all_nonself_pairs,
        )

        self.assertEqual(len(true_edges), 5)
        self.assertEqual(len(hub_candidates), 9)
        self.assertEqual(len(nonhub_candidates), 3)
        self.assertEqual(set(degree_table["Gene"]), genes)

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
