import math
import statistics
from pathlib import Path
from typing import Dict, List

import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.transforms import blended_transform_factory

from BLPlot._heatmap import (
    _DARK_COL,
    _draw_section,
    _setup_heatmap_axes,
)
from BLPlot.plotter import Plotter, iter_datasets_with_runs, random_classifier_baseline

plt.rcParams["font.size"] = 12


def _load_auprc_ratios(
    dataset_ids: List[str],
    dataset_paths: List[Path],
    gt_paths: List[Path],
) -> Dict[str, Dict[str, float]]:
    """
    Compute median AUPRC ratio per algorithm per dataset.

    For each dataset, reads AUPRC.csv and divides each per-run value by the
    random classifier baseline for that dataset, then takes the median across
    runs. Datasets with missing CSV or missing ground truth are skipped.

    Parameters
    ----------
    dataset_ids : list[str]
        Ordered dataset identifiers.
    dataset_paths : list[Path]
        Output directory for each dataset (same order as dataset_ids).
    gt_paths : list[Path]
        Ground truth CSV path for each dataset (same order as dataset_ids).

    Returns
    -------
    dict[str, dict[str, float]]
        algo -> {dataset_id -> median AUPRC ratio}, nan where unavailable.
    """
    result: Dict[str, Dict[str, float]] = {}

    for dataset_id, dataset_path, gt_path in zip(dataset_ids, dataset_paths, gt_paths):
        csv_path = dataset_path / 'AUPRC.csv'
        if not csv_path.exists():
            print(f"Warning: {csv_path} not found, skipping.")
            continue
        if not gt_path.exists():
            print(f"Warning: {gt_path} not found, skipping.")
            continue

        baseline = random_classifier_baseline(gt_path)
        if math.isnan(baseline) or baseline == 0.0:
            print(f"Warning: random baseline undefined for {gt_path}, skipping.")
            continue

        df = pd.read_csv(csv_path, index_col=0)
        for algo in df.index:
            vals = [v for v in df.loc[algo].tolist() if not math.isnan(v)]
            ratio = statistics.median([v / baseline for v in vals]) if vals else float('nan')
            result.setdefault(str(algo), {})[dataset_id] = ratio

    return result


def _load_spearman(
    dataset_ids: List[str],
    dataset_paths: List[Path],
) -> Dict[str, Dict[str, float]]:
    """
    Load median Spearman stability per algorithm per dataset.

    Reads Spearman.csv from each dataset's output directory. Datasets with
    missing CSV are skipped with a warning.

    Parameters
    ----------
    dataset_ids : list[str]
        Ordered dataset identifiers.
    dataset_paths : list[Path]
        Output directory for each dataset (same order as dataset_ids).

    Returns
    -------
    dict[str, dict[str, float]]
        algo -> {dataset_id -> MedianSpearman}, nan where unavailable.
    """
    result: Dict[str, Dict[str, float]] = {}

    for dataset_id, dataset_path in zip(dataset_ids, dataset_paths):
        csv_path = dataset_path / 'Spearman.csv'
        if not csv_path.exists():
            print(f"Warning: {csv_path} not found, skipping.")
            continue

        df = pd.read_csv(csv_path, index_col=0)
        for algo in df.index:
            val = df.loc[algo, 'MedianSpearman']
            result.setdefault(str(algo), {})[dataset_id] = float(val)

    return result


class PlotSummaryHeatmap(Plotter):
    """
    Replicates Figure 2 of Pratapa et al. 2020 (BEELINE).

    Produces a two-section heatmap: median AUPRC ratios (left) and median
    Spearman stability scores (right). Algorithms (rows) are sorted by
    decreasing median AUPRC ratio; datasets are columns. Each cell contains a
    rounded square sized and colored by its value. Values above the random
    predictor baseline (ratio >= 1) are drawn full-size with their raw value
    as white text. Alternating row backgrounds aid readability. Writes
    Summary.pdf to the output directory.
    """

    def __call__(self, config: dict, output_dir: Path, root: Path) -> None:
        """
        Generate the summary heatmap and write it to output_dir/Summary.pdf.

        Parameters
        ----------
        config : dict
            Parsed YAML configuration.
        output_dir : Path
            Directory where Summary.pdf is written.
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

        rows = list(iter_datasets_with_runs(config, root))
        if not rows:
            print("No datasets found for summary heatmap.")
            return

        dataset_ids   = [r[0] for r in rows]
        dataset_paths = [r[1] for r in rows]
        gt_paths      = [r[2] for r in rows]

        auprc_ratios = _load_auprc_ratios(dataset_ids, dataset_paths, gt_paths)
        spearman     = _load_spearman(dataset_ids, dataset_paths)

        all_algos = sorted(set(auprc_ratios) | set(spearman))
        if not all_algos:
            print("No data found for summary heatmap.")
            return

        # Sort algorithms by decreasing median-of-medians AUPRC ratio.
        def _median_auprc_ratio(algo: str) -> float:
            vals = [v for v in auprc_ratios.get(algo, {}).values() if not math.isnan(v)]
            return statistics.median(vals) if vals else 0.0

        sorted_algos = sorted(all_algos, key=_median_auprc_ratio, reverse=True)
        n_algos    = len(sorted_algos)
        n_datasets = len(dataset_ids)

        # Build raw value arrays (n_algos × n_datasets), sorted_algos[0] = best.
        auprc_arr = np.array([
            [auprc_ratios.get(a, {}).get(d, float('nan')) for d in dataset_ids]
            for a in sorted_algos
        ])
        spear_arr = np.array([
            [spearman.get(a, {}).get(d, float('nan')) for d in dataset_ids]
            for a in sorted_algos
        ])

        # --- Figure layout (mirrors old_heatmap.py) ---
        # Section 1 cols: x = 1 .. n_datasets
        # Gap:            x = n_datasets + 1
        # Section 2 cols: x = n_datasets + 2 .. n_datasets * 2 + 1
        total_cols = n_datasets * 2 + 1
        pad        = 2
        height     = 7
        asp_ratio  = (total_cols + pad) / (n_algos + pad)
        fig_size   = (height * asp_ratio + 0.5, height)

        # Gridspec: row 0 = heatmap (spans both columns), row 1 = palette legends.
        fig = plt.figure(figsize=(fig_size[0], fig_size[1] + 0.5))
        gs = fig.add_gridspec(2, 2, height_ratios=[n_algos + pad, 0.5], hspace=0.05)
        ax = fig.add_subplot(gs[0, :])

        _setup_heatmap_axes(ax, n_algos, sorted_algos, total_cols, pad)

        row_trans = blended_transform_factory(ax.transAxes, ax.transData)
        for row_idx in range(n_algos):
            bg = (0.9, 0.9, 0.9) if row_idx % 2 == 0 else (1.0, 1.0, 1.0)
            ax.add_artist(patches.Rectangle(
                (-0.25, n_algos - row_idx - 0.5),
                width=1.25, height=1,
                transform=row_trans,
                clip_on=False,
                edgecolor=(1, 1, 1), facecolor=bg,
            ))

        # Color palettes: viridis for AUPRC, cool (blue) for stability.
        auprc_palette = sns.color_palette("viridis", 11)
        spear_palette = sns.cubehelix_palette(11, reverse=True)

        _draw_section(
            ax, auprc_arr, n_algos, n_datasets,
            col_x_start=1,
            palette=auprc_palette,
            rand_cutoff=0.0,
            dataset_ids=dataset_ids,
            section_label='Median AUPRC ratios',
        )
        _draw_section(
            ax, spear_arr, n_algos, n_datasets,
            col_x_start=n_datasets + 2,
            palette=spear_palette,
            rand_cutoff=0.0,
            dataset_ids=dataset_ids,
            section_label='Median stability scores',
            switch_text=False,
        )

        # Color range legend: one palette bar per section at the bottom.
        for levl_idx, palette in enumerate([auprc_palette, spear_palette]):
            legend_ax = fig.add_subplot(gs[1, levl_idx])
            legend_ax.imshow(
                np.arange(len(palette)).reshape(1, len(palette)),
                cmap=mpl.colors.ListedColormap(list(palette)),
                interpolation='nearest', aspect='auto',
            )
            legend_ax.yaxis.set_ticks_position('none')
            legend_ax.xaxis.set_ticks_position('none')
            legend_ax.set_yticklabels([])
            legend_ax.set_xticks([0.5, len(palette) - 2])
            legend_ax.set_xticklabels(['Low/Poor', 'High/Good'], fontsize=12)

        out_path = output_dir / 'Summary.pdf'
        plt.savefig(out_path, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved summary heatmap to {out_path}")
