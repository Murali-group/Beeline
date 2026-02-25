from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.metrics import auc, precision_recall_curve

from BLPlot.plotter import (
    Plotter,
    get_algo_ids,
    iter_datasets_with_runs,
    load_dataset_metric,
    make_box_figure,
    random_classifier_baseline,
    write_pdf,
)


def _make_pr_curve_figure(
    run_path: Path,
    gt_path: Path,
    algos: List[str],
    dataset_id: str,
) -> 'plt.Figure | None':
    """
    Build a precision-recall curve figure for all algorithms in a single run.

    Each algorithm is drawn as a separate line. Algorithms with missing
    rankedEdges.csv or no true-positive predictions are skipped. Returns None
    if no algorithm produced a valid curve.

    Parameters
    ----------
    run_path : Path
        Output directory for the run (contains per-algorithm subdirectories).
    gt_path : Path
        Path to the ground truth edge list CSV (columns Gene1, Gene2).
    algos : list[str]
        Algorithm IDs to plot, drawn in sorted order.
    dataset_id : str
        Dataset name used in the plot title.

    Returns
    -------
    plt.Figure or None
        The created figure, or None if no valid curves could be drawn.
    """
    if not isinstance(run_path, Path):
        raise TypeError(f"run_path must be Path, got {type(run_path)}")
    if not isinstance(gt_path, Path):
        raise TypeError(f"gt_path must be Path, got {type(gt_path)}")

    if not gt_path.exists():
        print(f"Warning: ground truth not found at {gt_path}, skipping.")
        return None

    gt_df = pd.read_csv(gt_path, header=0)
    true_edges = set(zip(gt_df['Gene1'], gt_df['Gene2']))

    # Compute AUPRC for each algorithm first, then sort descending so that
    # the best-performing algorithm appears first in the legend and is drawn
    # with the first palette colour.
    curves = []
    for algo in algos:
        edges_path = run_path / algo / 'rankedEdges.csv'
        if not edges_path.exists():
            continue

        df = pd.read_csv(edges_path, sep='\t', header=0)
        predicted = df[df['Gene1'] != df['Gene2']].copy()
        if predicted.empty:
            continue

        labels = [
            1 if (g1, g2) in true_edges else 0
            for g1, g2 in zip(predicted['Gene1'], predicted['Gene2'])
        ]
        scores = predicted['EdgeWeight'].values

        # PR curve is undefined when no positive examples appear
        if sum(labels) == 0:
            continue

        precision, recall, _ = precision_recall_curve(labels, scores)
        score = auc(recall, precision)
        curves.append((score, algo, recall, precision))

    if not curves:
        return None

    # Sort by AUPRC descending
    curves.sort(key=lambda x: x[0], reverse=True)
    colors = sns.color_palette("Set1", n_colors=len(curves))

    fig, ax = plt.subplots(figsize=(7, 5))

    for (score, algo, recall, precision), color in zip(curves, colors):
        ax.plot(recall, precision, label=f'{algo} (AUPRC={score:.3f})', color=color)

    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title(f'Precision-Recall — {dataset_id}')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    plt.tight_layout()
    return fig


def _make_figure(
    dataset_id: str,
    dataset_path: Path,
    gt_path: Path,
    runs: list,
    algos: List[str],
) -> 'plt.Figure | None':
    """
    Produce the appropriate AUPRC figure for one dataset.

    Returns a PR curve figure for single-run datasets and a box plot for
    multi-run datasets.

    Parameters
    ----------
    dataset_id : str
        Dataset label used in the plot title.
    dataset_path : Path
        Output directory containing AUPRC.csv for multi-run datasets.
    gt_path : Path
        Ground truth CSV path (used for the random baseline and curve plots).
    runs : list
        Run dicts from the config for this dataset.
    algos : list[str]
        Algorithm IDs to include.

    Returns
    -------
    plt.Figure or None
    """
    if len(runs) == 1:
        run_path = dataset_path / runs[0]['run_id']
        return _make_pr_curve_figure(run_path, gt_path, algos, dataset_id)

    baseline = random_classifier_baseline(gt_path) if gt_path.exists() else None
    values = load_dataset_metric(dataset_path, 'AUPRC.csv')
    return make_box_figure(values, f'AUPRC — {dataset_id}', 'AUPRC', rand_value=baseline)


class PlotAUPRC(Plotter):
    """
    Plotter that produces one AUPRC graphic per dataset.

    For datasets with a single run a precision-recall curve is drawn (one line
    per algorithm). For datasets with multiple runs a box plot is drawn (one box
    per algorithm, distribution across runs). All pages are written to a single
    AUPRC.pdf in the output directory.
    """

    def __call__(self, config: dict, output_dir: Path, root: Path) -> None:
        """
        Generate per-dataset AUPRC graphics and write all pages to AUPRC.pdf.

        Parameters
        ----------
        config : dict
            Parsed YAML configuration.
        output_dir : Path
            Directory where AUPRC.pdf is written.
        root : Path
            Working directory from which config paths are resolved.

        Returns
        -------
        None
        """
        if not isinstance(config, dict):
            raise TypeError(f"config must be dict, got {type(config)}")
        if not isinstance(output_dir, Path):
            raise TypeError(f"output_dir must be Path, got {type(output_dir)}")
        if not isinstance(root, Path):
            raise TypeError(f"root must be Path, got {type(root)}")

        algos = get_algo_ids(config)
        write_pdf(
            output_dir / 'AUPRC.pdf',
            (
                _make_figure(did, dp, gtp, runs, algos)
                for did, dp, gtp, runs in iter_datasets_with_runs(config, root)
            ),
            'AUPRC',
        )
