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

from BLPlot.plotter import Plotter, iter_datasets_with_runs, random_classifier_baseline

plt.rcParams["font.size"] = 12

# Color used for cells whose value falls below the random-predictor cutoff.
_DARK_COL = '#2E4053'


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


def _flat_square(ax, cx: float, cy: float, col) -> None:
    """
    Draw a fixed-size flat square cell centered at (cx, cy).

    Matches the 'hm' shape in old_plotter.py: a borderless Rectangle with a
    slight horizontal gap (width = 0.9) so adjacent cells do not merge visually.

    Parameters
    ----------
    ax : matplotlib Axes
    cx, cy : float
        Center coordinates.
    col : color
        Fill and edge color (no visible border).
    """
    patch = patches.Rectangle(
        (cx - 0.5, cy - 0.5),
        width=0.9, height=1.0,
        facecolor=col, edgecolor=col,
    )
    ax.add_artist(patch)


def _draw_section(
    ax,
    data: np.ndarray,
    n_algos: int,
    n_datasets: int,
    col_x_start: int,
    palette: list,
    rand_cutoff: float,
    dataset_ids: List[str],
    section_label: str,
    switch_text: bool = False,
) -> None:
    """
    Draw one heatmap section onto ax using flat square patches.

    Cells whose raw value is below rand_cutoff are drawn in _DARK_COL. All
    other cells are colored by column-wise min-max normalization of the raw
    values. The raw value is annotated as text only in the column-maximum and
    column-minimum cells. Text color follows old_plotter.py: for switch_text=
    False (AUPRC) max-cell text is black and min-cell text is white; with
    switch_text=True (stability) the colors are swapped.

    Parameters
    ----------
    ax : matplotlib Axes
        Target axes.
    data : np.ndarray
        Shape (n_algos, n_datasets), raw values. Best algorithm at row 0
        (drawn at the top of the figure).
    n_algos : int
        Number of algorithm rows.
    n_datasets : int
        Number of dataset columns in this section.
    col_x_start : int
        X coordinate of the first column centre (1-indexed).
    palette : list
        Seaborn-style list of 11 colours (index 0 = low, 10 = high).
    rand_cutoff : float
        Values below this threshold are shown in _DARK_COL.
    dataset_ids : list[str]
        Column label for each dataset.
    section_label : str
        Header text displayed above the column labels.
    switch_text : bool
        When True, swap text colors so max-cell text is white and min-cell
        text is black (appropriate for dark-high palettes like cubehelix).
    """
    # Section header — centred above the column group
    mid_x = col_x_start + (n_datasets - 1) / 2
    ax.text(
        mid_x, n_algos + 2, section_label,
        fontsize=12, ha='center', va='center',
        bbox=dict(boxstyle='round', ec=(1, 1, 1), fc=(1, 1, 1)),
    )

    # Thick separator between section title and column labels
    x0 = col_x_start - 0.5
    x1 = col_x_start + n_datasets - 0.5
    ax.plot([x0, x1], [n_algos + 1.5, n_algos + 1.5],
            color='black', linewidth=2.5, solid_capstyle='butt')

    text_col_max = 'white' if switch_text else 'black'
    text_col_min = 'black' if switch_text else 'white'

    for col_idx in range(n_datasets):
        cx = col_x_start + col_idx

        # Column header — horizontal, centred above data rows
        ax.text(
            cx, n_algos + 1, dataset_ids[col_idx],
            fontsize=12, rotation=0, ha='center', va='center',
            bbox=dict(boxstyle='round', ec=(1, 1, 1, 0), fc=(1, 1, 1, 0)),
        )

        # Per-column min/max for normalization (all non-NaN values)
        col_vals = [data[r, col_idx] for r in range(n_algos)
                    if not math.isnan(data[r, col_idx])]
        col_min = min(col_vals) if col_vals else 0.0
        col_max = max(col_vals) if col_vals else 1.0

        for row_idx in range(n_algos):
            # Best algorithm at row 0 is drawn at the top (y = n_algos).
            cy  = n_algos - row_idx
            raw = data[row_idx, col_idx]

            if math.isnan(raw):
                ax.text(cx, cy, '-', fontsize=12,
                        ha='center', va='center', color=_DARK_COL)
                continue

            if raw < rand_cutoff:
                _flat_square(ax, cx, cy, _DARK_COL)
                continue

            # Column-wise min-max normalization for color selection
            if col_max > col_min:
                norm_val = (raw - col_min) / (col_max - col_min)
            else:
                norm_val = 0.5

            col = palette[min(int(np.floor(norm_val * 10)), 10)]
            _flat_square(ax, cx, cy, col)

            # Annotate the maximum and minimum value cells with the raw value.
            # Use one decimal place for values >= 1.0, two otherwise.
            fmt = f'{raw:.1f}' if raw >= 1.0 else f'{raw:.2f}'
            if norm_val >= 1.0:
                ax.text(cx, cy, fmt, fontsize=12,
                        ha='center', va='center', color=text_col_max,
                        bbox=dict(boxstyle='round', ec=(1,1,1,0), fc=(1,1,1,0)))
            elif norm_val <= 0.0:
                ax.text(cx, cy, fmt, fontsize=12,
                        ha='center', va='center', color=text_col_min,
                        bbox=dict(boxstyle='round', ec=(1,1,1,0), fc=(1,1,1,0)))


class PlotSummaryHeatmap(Plotter):
    """
    Replicates Figure 2 of Pratapa et al. 2020 (BEELINE).

    Produces a two-section heatmap: median AUPRC ratios (left) and median
    Spearman stability scores (right). Algorithms (rows) are sorted by
    decreasing median AUPRC ratio; datasets are columns. Each cell contains a
    rounded square sized and colored by its value. Values above the random
    predictor baseline (ratio >= 1) are drawn full-size with their raw value
    as white text. Alternating row backgrounds aid readability. Writes
    Overview.pdf to the output directory.
    """

    def __call__(self, config: dict, output_dir: Path, root: Path) -> None:
        """
        Generate the summary heatmap and write it to output_dir/Overview.pdf.

        Parameters
        ----------
        config : dict
            Parsed YAML configuration.
        output_dir : Path
            Directory where Overview.pdf is written.
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

        # Y-axis: algorithms listed top to bottom.
        # y = n_algos → top row (best), y = 1 → bottom row (worst).
        ax.set_yticks(np.arange(0, n_algos + pad))
        # Labels are listed bottom-to-top; best algo (sorted_algos[0]) is drawn
        # at the top (y = n_algos), so the label list must be reversed.
        ax.set_yticklabels([""] + sorted_algos[::-1] + [""], fontsize=12)

        # X-axis: no visible tick labels (headers drawn as text).
        ax.set_xticks(np.arange(0, total_cols + pad))
        ax.set_xticklabels([])
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')

        # Hide all spines.
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Alternating grey / white row backgrounds spanning both sections.
        for row_idx in range(n_algos):
            bg = (0.9, 0.9, 0.9) if row_idx % 2 == 0 else (1.0, 1.0, 1.0)
            ax.add_artist(patches.Rectangle(
                (0, n_algos - row_idx - 0.5),
                width=total_cols + 1, height=1,
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

        ax.set_xlim(0, total_cols + 1)
        ax.set_ylim(0, n_algos + 3)

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

        out_path = output_dir / 'Overview.pdf'
        plt.savefig(out_path, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved summary heatmap to {out_path}")
