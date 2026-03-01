#!/usr/bin/env python3
"""
plotDistributions.py – Plot per-gene pseudotime expression distributions
produced by ExperimentalRunner.

Each gene is drawn as a filled step curve over the [0, 1] pseudotime axis.
The y-axis shows the probability density (area under each curve = 1.0).

Usage:
    python BLPlot/plotDistributions.py DISTRIBUTIONS [-o OUTPUT]

Arguments:
    DISTRIBUTIONS   Path to distributions.txt (tab-separated; rows = genes,
                    columns = bin-centre pseudotime values).
    -o / --output   Output file path (PDF or PNG).
                    Default: distributions.pdf next to the input file.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def load_distributions(path: Path) -> pd.DataFrame:
    """
    Load a distributions.txt file produced by ExperimentalRunner.

    Parameters
    ----------
    path : Path
        Path to the tab-separated file. Rows are genes; columns are
        bin-centre pseudotime values (floats stored as strings).

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by gene name with float column names (bin centres)
        and float values (probability densities).

    Raises
    ------
    TypeError  – if path is not a Path.
    ValueError – if the file has fewer than 2 columns after the index.
    """
    if not isinstance(path, Path):
        raise TypeError(f"path must be Path, got {type(path)}")
    df = pd.read_csv(path, sep='\t', index_col=0)
    df.columns = df.columns.astype(float)
    if df.shape[1] < 2:
        raise ValueError(
            f"Expected at least 2 bin columns, got {df.shape[1]} in {path}")
    return df


def plot_distributions(df: pd.DataFrame, title: str = 'Gene Expression Distributions') -> plt.Figure:
    """
    Plot per-gene pseudotime expression distributions as one subplot per gene.

    Each gene gets its own axes with a filled step curve. Subplots are arranged
    in a grid and share the same x-axis range [0, 1] and y-axis scale so genes
    can be compared directly.

    Parameters
    ----------
    df : pd.DataFrame
        Distribution matrix; rows = genes, columns = bin-centre pseudotime values.
    title : str
        Figure suptitle.

    Returns
    -------
    plt.Figure
        The completed matplotlib figure.

    Raises
    ------
    TypeError  – if df is not a DataFrame or title is not a str.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"df must be DataFrame, got {type(df)}")
    if not isinstance(title, str):
        raise TypeError(f"title must be str, got {type(title)}")

    bin_centres = df.columns.values
    bin_width   = bin_centres[1] - bin_centres[0] if len(bin_centres) > 1 else 1.0
    left_edges  = bin_centres - bin_width / 2

    n_genes  = len(df)
    n_cols   = min(3, n_genes)
    n_rows   = (n_genes + n_cols - 1) // n_cols
    colors   = plt.rcParams['axes.prop_cycle'].by_key()['color']
    y_max    = df.values.max() * 1.05  # shared y ceiling across all subplots

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(5 * n_cols, 3 * n_rows),
                             sharey=True)
    # Flatten axes into a 1-D list regardless of grid shape.
    axes_flat = axes.flat if n_genes > 1 else [axes]

    for i, (gene, row) in enumerate(df.iterrows()):
        ax     = axes_flat[i]
        color  = colors[i % len(colors)]
        values = row.values

        ax.step(left_edges, values, where='post', color=color, linewidth=1.5)
        ax.fill_between(left_edges, values, step='post', color=color, alpha=0.25)

        ax.set_title(str(gene), fontsize=11)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, y_max)
        ax.set_xlabel('Pseudotime', fontsize=9)
        ax.set_ylabel('Prob. Density', fontsize=9)

    # Hide any unused subplot slots in the last row.
    for j in range(i + 1, n_rows * n_cols):
        axes_flat[j].set_visible(False)

    fig.suptitle(title, fontsize=13, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            'Plot per-gene pseudotime expression distributions from '
            'ExperimentalRunner output.'
        )
    )
    parser.add_argument(
        'distributions', type=Path,
        help='Path to distributions.txt (tab-separated: rows=genes, cols=bin centres)',
    )
    parser.add_argument(
        '-o', '--output', type=Path, default=None,
        help='Output file path (PDF or PNG). Default: distributions.pdf beside the input.',
    )
    args = parser.parse_args()

    if not args.distributions.exists():
        print(f"Error: file not found: {args.distributions}", file=sys.stderr)
        sys.exit(1)

    output = args.output or args.distributions.parent / 'distributions.pdf'

    df    = load_distributions(args.distributions)
    title = f'Gene Expression Distributions — {args.distributions.parent.parent.name}'
    fig   = plot_distributions(df, title=title)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {output}")


if __name__ == '__main__':
    main()
