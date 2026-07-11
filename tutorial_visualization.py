#!/usr/bin/env python
"""
tutorial_visualization.py

Draws a single GRN inference result as an undirected network graph so tutorial
participants can visually inspect how one algorithm's predictions compare to a
ground truth network.

The top-n edges of the ranked edge list are taken, where n is the number of
edges in the ground truth network. Each drawn edge is coloured by category:
  * green  - correct prediction   (in top-n AND in ground truth)
  * red    - incorrect prediction (in top-n but NOT in ground truth)
  * grey   - unpredicted GT edge   (in ground truth but NOT in top-n)
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")  # headless backend so the script works over SSH / in CI
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# Colours for each edge category.
COLOR_CORRECT = "#2ca02c"   # green
COLOR_INCORRECT = "#d62728"  # red
COLOR_MISSED = "#b0b0b0"    # grey


def load_edges(path, weighted):
    """
    Read an edge list into an ordered list of (Gene1, Gene2) tuples.

    Parameters:
        path (str or Path): Path to a CSV/TSV edge list whose first two columns
            are the source and target gene. The delimiter (comma or tab) is
            sniffed automatically.
        weighted (bool): If True the third column is treated as an edge weight
            and rows are returned sorted by descending weight (highest-ranked
            first). If False rows are returned in file order.

    Returns:
        list of (str, str): Directed edges, ordered as described above.
    """
    if not Path(path).is_file():
        sys.exit("Error: file not found: {}".format(path))

    # sep=None with the python engine sniffs comma vs. tab delimiters.
    df = pd.read_csv(path, sep=None, engine="python")

    if df.shape[1] < 2:
        sys.exit("Error: {} must have at least two columns (Gene1, Gene2)".format(path))

    if weighted and df.shape[1] >= 3:
        weight_col = df.columns[2]
        df = df.sort_values(by=weight_col, ascending=False)

    g1, g2 = df.columns[0], df.columns[1]
    return [(str(a), str(b)) for a, b in zip(df[g1], df[g2])]


def main():
    parser = argparse.ArgumentParser(
        description="Visualize a single GRN inference result against a ground truth network."
    )
    parser.add_argument("ground_truth", help="Path to the ground truth network CSV (Gene1, Gene2).")
    parser.add_argument("ranked_edges", help="Path to the ranked edge list (Gene1, Gene2, EdgeWeight).")
    parser.add_argument(
        "-o", "--output", default="tutorial_visualization.png",
        help="Output image path (default: tutorial_visualization.png).",
    )
    parser.add_argument(
        "-n", "--top-n", type=int, default=None,
        help="Override the number of top ranked edges to draw "
             "(default: number of edges in the ground truth).",
    )
    parser.add_argument(
        "--seed", type=int, default=0,
        help="Random seed for the layout, for reproducible plots (default: 0).",
    )
    args = parser.parse_args()

    # Ground truth is unweighted; ranked edges are weighted and must be sorted.
    gt_edges = load_edges(args.ground_truth, weighted=False)
    ranked_edges = load_edges(args.ranked_edges, weighted=True)

    # Edges are undirected: canonicalise each as a sorted (Gene1, Gene2) tuple
    # so that A-B and B-A collapse to the same edge.
    def canon(edge):
        a, b = edge
        return (a, b) if a <= b else (b, a)

    gt_set = {canon(e) for e in gt_edges}
    n = args.top_n if args.top_n is not None else len(gt_set)

    # Walk the ranked list in order, keeping each undirected edge only once
    # (its highest-ranked occurrence), then take the top-n.
    seen = set()
    pred_edges = []
    for e in ranked_edges:
        c = canon(e)
        if c not in seen:
            seen.add(c)
            pred_edges.append(c)
    pred_edges = pred_edges[:n]
    pred_set = set(pred_edges)

    # Classify edges into the three categories.
    correct = pred_set & gt_set          # green: predicted and true
    incorrect = pred_set - gt_set        # red:   predicted but not true
    missed = gt_set - pred_set           # grey:  true but not predicted

    print("Ground truth edges : {}".format(len(gt_set)))
    print("Top-n predictions  : {}".format(len(pred_set)))
    print("  correct (green)  : {}".format(len(correct)))
    print("  incorrect (red)  : {}".format(len(incorrect)))
    print("  missed (grey)    : {}".format(len(missed)))

    # Build a directed graph containing every edge we will draw, keeping the
    # three categories in separate lists so they can be coloured independently.
    graph = nx.Graph()
    graph.add_edges_from(correct)
    graph.add_edges_from(incorrect)
    graph.add_edges_from(missed)

    # Scale node/label size to graph size so large networks stay legible.
    n_nodes = graph.number_of_nodes()
    node_size = max(130, min(5200, 104000 / max(n_nodes, 1)))
    # Scale label font with node size (larger nodes -> larger labels).
    label_font = max(8, min(18, int(node_size ** 0.5 / 3)))
    with_labels = n_nodes <= 60

    # Arrange nodes evenly around a ring. sorted() gives a stable, readable order.
    pos = nx.circular_layout(graph.subgraph(sorted(graph.nodes())))

    fig, ax = plt.subplots(figsize=(14, 14))

    # Draw edges FIRST, then nodes on top, so edges are cleanly hidden behind the
    # node discs at their endpoints instead of clipping across them.
    # Grey (missed) is drawn first so correct/incorrect edges sit visually on top.
    for edge_set, color in ((missed, COLOR_MISSED),
                            (incorrect, COLOR_INCORRECT),
                            (correct, COLOR_CORRECT)):
        if edge_set:
            nx.draw_networkx_edges(
                graph, pos, edgelist=list(edge_set), edge_color=color,
                width=1.2, arrows=False, ax=ax,
            )

    nx.draw_networkx_nodes(graph, pos, node_size=node_size,
                           node_color="#cfd8e3", edgecolors="#555555",
                           linewidths=0.8, ax=ax)

    if with_labels:
        nx.draw_networkx_labels(graph, pos, font_size=label_font, ax=ax)

    legend_handles = [
        Line2D([0], [0], color=COLOR_CORRECT, lw=2,
               label="Correct prediction ({})".format(len(correct))),
        Line2D([0], [0], color=COLOR_INCORRECT, lw=2,
               label="Incorrect prediction ({})".format(len(incorrect))),
        Line2D([0], [0], color=COLOR_MISSED, lw=2,
               label="Unpredicted GT edge ({})".format(len(missed))),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=12)

    ax.set_title(
        "GRN prediction vs. ground truth (top {} edges)".format(n), fontsize=16
    )
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150, bbox_inches="tight")
    print("Saved figure to {}".format(args.output))


if __name__ == "__main__":
    main()
