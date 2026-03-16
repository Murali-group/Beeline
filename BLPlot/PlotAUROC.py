from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.metrics import auc, roc_curve

from BLPlot.plotter import (
    Plotter,
    get_algo_ids,
    iter_datasets_with_runs,
    load_dataset_metric,
    make_box_figure,
)


def _make_roc_curve_figure(
    run_path: Path,
    gt_path: Path,
    algos: List[str],
    dataset_label: str,
) -> 'plt.Figure | None':
    """
    Build a ROC curve figure for all algorithms in a single run.

    Each algorithm is drawn as a separate line. Algorithms with missing
    rankedEdges.csv or only one class in predictions are skipped. Returns
    None if no algorithm produced a valid curve.

    Parameters
    ----------
    run_path : Path
        Output directory for the run (contains per-algorithm subdirectories).
    gt_path : Path
        Path to the ground truth edge list CSV (columns Gene1, Gene2).
    algos : list[str]
        Algorithm IDs to plot, drawn in sorted order.
    dataset_label : str
        Dataset label used in the plot title (nickname if set, else dataset_id).

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

    # Compute AUROC for each algorithm first, then sort descending so that
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

        # ROC curve is undefined when only one class appears
        if sum(labels) == 0 or sum(labels) == len(labels):
            continue

        fpr, tpr, _ = roc_curve(labels, scores)
        score = auc(fpr, tpr)
        curves.append((score, algo, fpr, tpr))

    if not curves:
        return None

    # Sort by AUROC descending
    curves.sort(key=lambda x: x[0], reverse=True)
    colors = sns.color_palette("Set1", n_colors=len(curves))

    fig, ax = plt.subplots(figsize=(7, 5))
    # Random classifier diagonal reference line
    ax.plot([0, 1], [0, 1], color='grey', linestyle='--', linewidth=0.8, label='Random')

    for (score, algo, fpr, tpr), color in zip(curves, colors):
        ax.plot(fpr, tpr, label=f'{algo} (AUROC={score:.3f})', color=color)

    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(f'ROC Curve — {dataset_label}')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    plt.tight_layout()
    return fig


def _make_figure(
    dataset_label: str,
    dataset_path: Path,
    gt_path: Path,
    runs: list,
    algos: List[str],
) -> 'plt.Figure | None':
    """
    Produce the appropriate AUROC figure for one dataset.

    Returns a ROC curve figure for single-run datasets and a box plot for
    multi-run datasets. The random classifier AUROC baseline is always 0.5.

    Parameters
    ----------
    dataset_label : str
        Dataset label used in the plot title (nickname if set, else dataset_id).
    dataset_path : Path
        Output directory containing AUROC.csv for multi-run datasets.
    gt_path : Path
        Ground truth CSV path (used for single-run curve plots).
    runs : list
        Run dicts from the config for this dataset.
    algos : list[str]
        Algorithm IDs to include.

    Returns
    -------
    plt.Figure or None
    """
    if len(runs) == 1:
        # In single_run mode run_id is None; output lives directly under dataset_path.
        run_id = runs[0]['run_id']
        run_path = dataset_path / run_id if run_id is not None else dataset_path
        return _make_roc_curve_figure(run_path, gt_path, algos, dataset_label)

    values = load_dataset_metric(dataset_path, 'AUROC.csv')
    return make_box_figure(values, f'AUROC — {dataset_label}', 'AUROC', rand_value=0.5)


class PlotAUROC(Plotter):
    """
    Plotter that produces one AUROC graphic per dataset.

    For datasets with a single run a ROC curve is drawn (one line per
    algorithm). For datasets with multiple runs a box plot is drawn (one box
    per algorithm, distribution across runs). Each dataset is written as both a
    PDF and PNG under an AUROC/ subdirectory of the output directory, named
    <dataset_label>-AUROC.pdf and <dataset_label>-AUROC.png.
    """

    def __call__(self, config: dict, output_dir: Path, root: Path) -> None:
        """
        Generate per-dataset AUROC graphics, writing each to its own PDF and PNG.

        Creates output_dir/AUROC/ and writes <dataset_label>-AUROC.pdf and
        <dataset_label>-AUROC.png for each enabled dataset. Datasets that
        produce no figure (missing data or no valid curves) are skipped.

        Parameters
        ----------
        config : dict
            Parsed YAML configuration.
        output_dir : Path
            Parent directory; an AUROC/ subdirectory is created inside it.
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

        auroc_dir = output_dir / 'AUROC'
        auroc_dir.mkdir(parents=True, exist_ok=True)

        algos = get_algo_ids(config)
        for _, dlabel, dp, gtp, runs in iter_datasets_with_runs(config, root):
            fig = _make_figure(dlabel, dp, gtp, runs, algos)
            if fig is None:
                continue
            safe_label = dlabel.replace('/', '-')
            stem = auroc_dir / f'{safe_label}-AUROC'
            fig.savefig(stem.with_suffix('.pdf'))
            fig.savefig(stem.with_suffix('.png'))
            plt.close(fig)
            print(f"Saved {stem}.pdf and .png")
