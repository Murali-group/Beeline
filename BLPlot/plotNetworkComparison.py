#!/usr/bin/env python3
"""
plotNetworkComparison.py – Plot one or more ranked-edges networks beside a
ground-truth network, using identical node positions so the graphs can be
compared directly.

Usage:
    python BLPlot/plotNetworkComparison.py RANKED_EDGES [RANKED_EDGES ...] -g GROUND_TRUTH -k K [-o OUTPUT]

Arguments:
    RANKED_EDGES         One or more paths to rankedEdges.csv files
                         (tab-separated: Gene1, Gene2, EdgeWeight).
                         Each file produces one panel. Panel titles are taken
                         from the immediate parent directory of each file
                         (e.g. the algorithm name in BEELINE's output structure).
    -g / --ground-truth  Path to GroundTruthNetwork.csv (comma-separated: Gene1, Gene2, Type)
    -k / --top-k         Number of top-ranked edges to show per ranked-edges panel.
    -o / --output        Output file path (PDF or PNG). Default: networkComparison.pdf
"""

import argparse
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


# -- Colours ------------------------------------------------------------------

# Ground-truth edge colours
_ACTIVATE_COLOR = '#2ca02c'     # green  – activation (+)
_REPRESS_COLOR  = '#d62728'     # red    – repression (-)
_OTHER_COLOR    = '#7f7f7f'     # grey   – any other sign value

# Ranked-edges panel node colours
_CONNECTED_NODE_COLOR  = '#aec7e8'  # light blue – gene with at least one edge
_ISOLATED_NODE_COLOR   = '#d9d9d9'  # grey  – gene present but unconnected
# Ranked-edges use the same green/red as the ground-truth panel for edge sign.


# -- Data loading -------------------------------------------------------------

def load_ground_truth(path: Path) -> pd.DataFrame:
    """
    Load a GroundTruthNetwork.csv file.

    Parameters
    ----------
    path : Path
        Path to a comma-separated file with columns Gene1, Gene2, Type.

    Returns
    -------
    pd.DataFrame
        DataFrame with string columns Gene1, Gene2, Type.

    Raises
    ------
    TypeError  – if path is not a Path object.
    ValueError – if required columns are missing.
    """
    if not isinstance(path, Path):
        raise TypeError(f"path must be Path, got {type(path)}")
    df = pd.read_csv(path)
    required = {'Gene1', 'Gene2', 'Type'}
    if not required.issubset(df.columns):
        raise ValueError(
            f"Ground truth CSV must have columns {required}, got {list(df.columns)}"
        )
    df['Gene1'] = df['Gene1'].astype(str)
    df['Gene2'] = df['Gene2'].astype(str)
    df['Type']  = df['Type'].astype(str)
    return df


def load_ranked_edges(path: Path, top_k: int) -> pd.DataFrame:
    """
    Load a rankedEdges.csv file and return the top-K rows by EdgeWeight.

    Parameters
    ----------
    path : Path
        Path to a tab-separated file with columns Gene1, Gene2, EdgeWeight.
    top_k : int
        Maximum number of edges to return (highest-weight first).

    Returns
    -------
    pd.DataFrame
        DataFrame with string columns Gene1, Gene2 and float column EdgeWeight,
        sorted descending by absolute EdgeWeight, length at most top_k.

    Raises
    ------
    TypeError  – if path is not a Path or top_k is not an int.
    ValueError – if required columns are missing or top_k < 1.
    """
    if not isinstance(path, Path):
        raise TypeError(f"path must be Path, got {type(path)}")
    if not isinstance(top_k, int) or top_k < 1:
        raise ValueError(f"top_k must be a positive int, got {top_k!r}")
    df = pd.read_csv(path, sep='\t')
    required = {'Gene1', 'Gene2', 'EdgeWeight'}
    if not required.issubset(df.columns):
        raise ValueError(
            f"Ranked edges CSV must have columns {required}, got {list(df.columns)}"
        )
    df['Gene1']      = df['Gene1'].astype(str)
    df['Gene2']      = df['Gene2'].astype(str)
    df['EdgeWeight'] = df['EdgeWeight'].astype(float)
    # Remove self-edges before ranking so top_k counts only real edges.
    df = df[df['Gene1'] != df['Gene2']]
    df = df.reindex(df['EdgeWeight'].abs().sort_values(ascending=False).index).head(top_k)
    return df


# -- Graph construction -------------------------------------------------------

def _build_gt_graph(gt_df: pd.DataFrame) -> nx.DiGraph:
    """
    Build a directed graph from a ground-truth DataFrame.

    Parameters
    ----------
    gt_df : pd.DataFrame
        DataFrame with columns Gene1, Gene2, Type.

    Returns
    -------
    nx.DiGraph
        Directed graph; each edge carries a 'sign' attribute ('+', '-', or other).
    """
    if not isinstance(gt_df, pd.DataFrame):
        raise TypeError(f"gt_df must be DataFrame, got {type(gt_df)}")
    G = nx.DiGraph()
    for _, row in gt_df.iterrows():
        if row['Gene1'] != row['Gene2']:
            G.add_edge(row['Gene1'], row['Gene2'], sign=row['Type'])
    return G


def _build_ranked_graph(ranked_df: pd.DataFrame, nodes: set) -> nx.DiGraph:
    """
    Build a directed graph from a ranked-edges DataFrame.

    All nodes in `nodes` are added to the graph before edges are inserted,
    so genes from the ground truth that have no ranked edges are still present
    as isolated nodes.

    Parameters
    ----------
    ranked_df : pd.DataFrame
        DataFrame with columns Gene1, Gene2, EdgeWeight.
    nodes : set
        Full set of gene names to include in the graph.

    Returns
    -------
    nx.DiGraph
        Directed graph; each edge carries a 'weight' attribute (float).
    """
    if not isinstance(ranked_df, pd.DataFrame):
        raise TypeError(f"ranked_df must be DataFrame, got {type(ranked_df)}")
    if not isinstance(nodes, set):
        raise TypeError(f"nodes must be set, got {type(nodes)}")
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    for _, row in ranked_df.iterrows():
        if row['Gene1'] != row['Gene2']:
            G.add_edge(row['Gene1'], row['Gene2'], weight=float(row['EdgeWeight']))
    return G


# -- Drawing helpers ----------------------------------------------------------

# Node radius used for arrow shrink calculations.  Matches node_size=600
# passed to draw_networkx_nodes (matplotlib scatter area in pt^2 → radius ≈ √(600/π) ≈ 14 pt).
# A small extra margin is added so arrowheads do not clip into the node circle.
_NODE_RADIUS_PT = 16


def _draw_edges(ax: plt.Axes, edgelist: list, pos: dict,
                color: str, lw: float, rad: float = 0.2) -> None:
    """
    Draw directed edges on `ax` using matplotlib annotate arrows.

    Each edge in `edgelist` is drawn as a FancyArrowPatch via ax.annotate.
    The `rad` parameter controls the curvature of the arc; a non-zero value is
    used to separate both members of a bidirectional pair so they are visually
    distinct.  Unlike networkx's draw_networkx_edges, this approach supports
    connectionstyle in all networkx versions.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    edgelist : list of (node, node)
        Edges to draw.
    pos : dict
        Node positions (data coordinates) from _compute_layout.
    color : str
        Edge and arrowhead colour.
    lw : float
        Line width in points.
    rad : float
        Arc radius passed to connectionstyle 'arc3,rad=<rad>'.  Zero gives a
        straight arrow; positive curves counter-clockwise (from source to target).

    Returns
    -------
    None
    """
    if not isinstance(ax, plt.Axes):
        raise TypeError(f"ax must be plt.Axes, got {type(ax)}")
    if not isinstance(pos, dict):
        raise TypeError(f"pos must be dict, got {type(pos)}")

    for u, v in edgelist:
        ax.annotate(
            "",
            xy=pos[v], xycoords='data',
            xytext=pos[u], textcoords='data',
            arrowprops=dict(
                arrowstyle='->',
                connectionstyle=f'arc3,rad={rad}',
                color=color,
                lw=lw,
                # Shrink arrows back from node centres so heads stop at the
                # node boundary rather than overlapping the node circle.
                shrinkA=_NODE_RADIUS_PT,
                shrinkB=_NODE_RADIUS_PT,
            ),
        )


def _compute_layout(G: nx.DiGraph) -> dict:
    """
    Compute a Kamada-Kawai layout for the graph.

    Kamada-Kawai minimises graph-theoretic distances between nodes and
    produces evenly spread, non-collinear arrangements more reliably than
    spring layout for small graphs.

    Parameters
    ----------
    G : nx.DiGraph
        Graph whose nodes define the layout.

    Returns
    -------
    dict
        Mapping from node name to (x, y) numpy array.
    """
    if not isinstance(G, nx.DiGraph):
        raise TypeError(f"G must be nx.DiGraph, got {type(G)}")
    return nx.kamada_kawai_layout(G)


def _draw_gt_panel(ax: plt.Axes, G: nx.DiGraph, pos: dict) -> None:
    """
    Draw the ground-truth graph on `ax`, colouring edges by sign.

    Activation edges (+) are drawn green; repression edges (-) are drawn red;
    all other sign values are drawn grey. Nodes are drawn in light blue.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    G : nx.DiGraph
        Ground-truth directed graph with a 'sign' edge attribute.
    pos : dict
        Node positions returned by _compute_layout.

    Returns
    -------
    None
    """
    if not isinstance(ax, plt.Axes):
        raise TypeError(f"ax must be plt.Axes, got {type(ax)}")
    if not isinstance(G, nx.DiGraph):
        raise TypeError(f"G must be nx.DiGraph, got {type(G)}")
    if not isinstance(pos, dict):
        raise TypeError(f"pos must be dict, got {type(pos)}")

    activate = [(u, v) for u, v, d in G.edges(data=True) if d.get('sign') == '+']
    repress  = [(u, v) for u, v, d in G.edges(data=True) if d.get('sign') == '-']
    other    = [(u, v) for u, v, d in G.edges(data=True)
                if d.get('sign') not in ('+', '-')]

    nx.draw_networkx_nodes(G, pos, ax=ax,
                           node_color=_CONNECTED_NODE_COLOR, node_size=600)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=9)

    # Curve both members of a bidirectional pair so they are visually distinct.
    edge_set = set(G.edges())
    def _rad(u, v):
        return 0.2 if (v, u) in edge_set else 0.0

    for u, v in activate:
        _draw_edges(ax, [(u, v)], pos, _ACTIVATE_COLOR, lw=2.0, rad=_rad(u, v))
    for u, v in repress:
        _draw_edges(ax, [(u, v)], pos, _REPRESS_COLOR,  lw=2.0, rad=_rad(u, v))
    for u, v in other:
        _draw_edges(ax, [(u, v)], pos, _OTHER_COLOR,    lw=2.0, rad=_rad(u, v))

    ax.legend(handles=[
        mpatches.Patch(color=_ACTIVATE_COLOR, label='Activation (+)'),
        mpatches.Patch(color=_REPRESS_COLOR,  label='Repression (−)'),
    ], loc='upper left', fontsize=8)


def _draw_ranked_panel(ax: plt.Axes, G: nx.DiGraph, pos: dict,
                       top_k: int) -> None:
    """
    Draw the ranked-edges graph on `ax`.

    Edges with a positive weight are drawn in green (matching the ground-truth
    activation colour); edges with a negative weight are drawn in red (matching
    repression); edges with no weight attribute default to grey.  Connected nodes
    are drawn in light blue; isolated nodes are drawn in grey.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    G : nx.DiGraph
        Ranked-edges directed graph with all ground-truth genes pre-populated.
        Each edge should carry a 'weight' attribute (float).
    pos : dict
        Node positions returned by _compute_layout.
    top_k : int
        Actual number of edges shown (used for the panel title).

    Returns
    -------
    None
    """
    if not isinstance(ax, plt.Axes):
        raise TypeError(f"ax must be plt.Axes, got {type(ax)}")
    if not isinstance(G, nx.DiGraph):
        raise TypeError(f"G must be nx.DiGraph, got {type(G)}")
    if not isinstance(pos, dict):
        raise TypeError(f"pos must be dict, got {type(pos)}")
    if not isinstance(top_k, int) or top_k < 0:
        raise ValueError(f"top_k must be a non-negative int, got {top_k!r}")

    connected = {n for n in G.nodes() if G.degree(n) > 0}
    node_colors = [
        _CONNECTED_NODE_COLOR if n in connected else _ISOLATED_NODE_COLOR
        for n in G.nodes()
    ]

    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors, node_size=600)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=9)

    if G.edges():
        edge_set = set(G.edges())
        for u, v, data in G.edges(data=True):
            weight = data.get('weight', 0.0)
            if weight > 0:
                color = _ACTIVATE_COLOR
            elif weight < 0:
                color = _REPRESS_COLOR
            else:
                color = _OTHER_COLOR
            rad = 0.1 if (v, u) in edge_set else 0.0
            _draw_edges(ax, [(u, v)], pos, color, lw=1.5, rad=rad)

    ax.legend(handles=[
        mpatches.Patch(color=_ACTIVATE_COLOR,      label='Positive weight'),
        mpatches.Patch(color=_REPRESS_COLOR,        label='Negative weight'),
        mpatches.Patch(color=_CONNECTED_NODE_COLOR, label='Connected gene'),
        mpatches.Patch(color=_ISOLATED_NODE_COLOR,  label='Isolated gene'),
    ], loc='upper left', fontsize=8)


# -- Figure assembly ----------------------------------------------------------

def make_comparison_figure(gt_df: pd.DataFrame,
                           ranked_panels: list) -> plt.Figure:
    """
    Build the network comparison figure with one ground-truth panel followed
    by one panel per ranked-edges input.

    Layout is computed from the ground-truth graph so its edge structure drives
    node spacing. All panels share identical node positions.

    Parameters
    ----------
    gt_df : pd.DataFrame
        Ground-truth DataFrame with columns Gene1, Gene2, Type.
    ranked_panels : list of (pd.DataFrame, str, int)
        Each entry is a tuple of:
          - ranked_df  : ranked-edges DataFrame (Gene1, Gene2, EdgeWeight),
                         already trimmed to at most top_k rows.
          - label      : panel title string (e.g. algorithm name).
          - actual_k   : number of edges in ranked_df (may be less than
                         the requested top_k if the file had fewer rows).

    Returns
    -------
    plt.Figure
        The completed matplotlib figure.

    Raises
    ------
    TypeError  – if gt_df is not a DataFrame or ranked_panels is not a list.
    ValueError – if ranked_panels is empty.
    """
    if not isinstance(gt_df, pd.DataFrame):
        raise TypeError(f"gt_df must be DataFrame, got {type(gt_df)}")
    if not isinstance(ranked_panels, list):
        raise TypeError(f"ranked_panels must be list, got {type(ranked_panels)}")
    if not ranked_panels:
        raise ValueError("ranked_panels must contain at least one entry")

    gt_graph = _build_gt_graph(gt_df)
    gt_nodes = set(gt_graph.nodes())

    # Layout derived from the ground-truth graph so its edges drive spacing.
    pos = _compute_layout(gt_graph)

    n_panels = 1 + len(ranked_panels)
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 6))
    # Ensure axes is always a list even when n_panels == 1.
    if n_panels == 1:
        axes = [axes]
    fig.suptitle('Network Comparison', fontsize=14, fontweight='bold')

    ax_gt = axes[0]
    ax_gt.set_title('Ground Truth Network')
    _draw_gt_panel(ax_gt, gt_graph, pos)
    ax_gt.axis('off')

    for ax, (ranked_df, label, actual_k) in zip(axes[1:], ranked_panels):
        ranked_graph = _build_ranked_graph(ranked_df, gt_nodes)
        ax.set_title(f'{label} — top {actual_k}')
        _draw_ranked_panel(ax, ranked_graph, pos, actual_k)
        ax.axis('off')

    # rect=[left, bottom, right, top] reserves the top 8% for the suptitle
    # so it does not overlap with panel titles.
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    return fig


# -- Entry point --------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            'Plot a ranked-edges network beside a ground-truth network with '
            'identical node positions. All genes from the ground truth appear '
            'in the ranked-edges panel even when isolated.'
        )
    )
    parser.add_argument(
        'ranked_edges', type=Path, nargs='+',
        help=(
            'One or more paths to rankedEdges.csv files '
            '(tab-separated: Gene1, Gene2, EdgeWeight). '
            'Each produces one panel; the panel title is taken from the '
            "file's immediate parent directory name."
        ),
    )
    parser.add_argument(
        '-g', '--ground-truth', type=Path, required=True,
        dest='ground_truth',
        help='Path to GroundTruthNetwork.csv (comma-separated: Gene1, Gene2, Type)',
    )
    parser.add_argument(
        '-k', '--top-k', type=int, required=True,
        help='Number of top-ranked edges to display per panel',
    )
    parser.add_argument(
        '-o', '--output', type=Path, default=Path('networkComparison.pdf'),
        help='Output file path (PDF or PNG). Default: networkComparison.pdf',
    )
    args = parser.parse_args()

    errors = False
    for p in args.ranked_edges:
        if not p.exists():
            print(f"Error: ranked edges file not found: {p}", file=sys.stderr)
            errors = True
    if not args.ground_truth.exists():
        print(f"Error: ground truth file not found: {args.ground_truth}",
              file=sys.stderr)
        errors = True
    if args.top_k < 1:
        print(f"Error: --top-k must be at least 1, got {args.top_k}",
              file=sys.stderr)
        errors = True
    if errors:
        sys.exit(1)

    gt_df = load_ground_truth(args.ground_truth)

    ranked_panels = []
    for p in args.ranked_edges:
        # Use the immediate parent directory as the label (e.g. algorithm name).
        label     = p.parent.name
        ranked_df = load_ranked_edges(p, args.top_k)
        actual_k  = len(ranked_df)
        if actual_k < args.top_k:
            print(f"Warning: {label}: only {actual_k} edges available "
                  f"(requested {args.top_k})")
        ranked_panels.append((ranked_df, label, actual_k))

    fig = make_comparison_figure(gt_df, ranked_panels)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {args.output}")


if __name__ == '__main__':
    main()
