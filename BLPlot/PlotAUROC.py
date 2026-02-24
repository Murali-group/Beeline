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
    write_pdf,
)


def _make_roc_curve_figure(
    run_path: Path,
    gt_path: Path,
    algos: List[str],
    dataset_id: str,
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

    sorted_algos = sorted(algos)
    colors = sns.color_palette("Set1", n_colors=len(sorted_algos))

    fig, ax = plt.subplots(figsize=(7, 5))
    # Random classifier diagonal reference line
    ax.plot([0, 1], [0, 1], color='grey', linestyle='--', linewidth=0.8, label='Random')
    any_line = False

    for algo, color in zip(sorted_algos, colors):
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
        ax.plot(fpr, tpr, label=f'{algo} (AUROC={score:.3f})', color=color)
        any_line = True

    if not any_line:
        plt.close(fig)
        return None

    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(f'ROC Curve — {dataset_id}')
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
    Produce the appropriate AUROC figure for one dataset.

    Returns a ROC curve figure for single-run datasets and a box plot for
    multi-run datasets. The random classifier AUROC baseline is always 0.5.

    Parameters
    ----------
    dataset_id : str
        Dataset label used in the plot title.
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
        run_path = dataset_path / runs[0]['run_id']
        return _make_roc_curve_figure(run_path, gt_path, algos, dataset_id)

    values = load_dataset_metric(dataset_path, 'AUROC.csv')
    return make_box_figure(values, f'AUROC — {dataset_id}', 'AUROC', rand_value=0.5)


class PlotAUROC(Plotter):
    """
    Plotter that produces one AUROC graphic per dataset.

    For datasets with a single run a ROC curve is drawn (one line per
    algorithm). For datasets with multiple runs a box plot is drawn (one box
    per algorithm, distribution across runs). All pages are written to a single
    AUROC.pdf in the output directory.
    """

    def __call__(self, config: dict, output_dir: Path, root: Path) -> None:
        """
        Generate per-dataset AUROC graphics and write all pages to AUROC.pdf.

        Parameters
        ----------
        config : dict
            Parsed YAML configuration.
        output_dir : Path
            Directory where AUROC.pdf is written.
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
            output_dir / 'AUROC.pdf',
            (
                _make_figure(did, dp, gtp, runs, algos)
                for did, dp, gtp, runs in iter_datasets_with_runs(config, root)
            ),
            'AUROC',
        )
