import math
import statistics
from pathlib import Path
from typing import Dict, List, Tuple

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
from BLPlot.plotter import Plotter, iter_datasets_with_runs

plt.rcParams["font.size"] = 12


def _signed_random_baselines(
    gt_path: Path,
) -> Tuple[float, float, float]:
    """
    Compute random-predictor EPR baselines for a dataset.

    For the full network and each signed edge type the baseline is
    k / (n*(n-1)), where k is the edge count for that type and n is the
    number of unique genes across the whole network.

    Parameters
    ----------
    gt_path : Path
        Ground truth CSV with columns Gene1, Gene2, and optionally Type
        ('+' for activation, '-' for inhibitory).

    Returns
    -------
    tuple of (baseline_total, baseline_activation, baseline_inhibitory)
        Random EPR baselines; a component is nan when its k is 0 or
        n_possible is 0.
    """
    if not isinstance(gt_path, Path):
        raise TypeError(f"gt_path must be Path, got {type(gt_path)}")

    gt = pd.read_csv(gt_path, header=0)
    gt = gt[gt['Gene1'] != gt['Gene2']]

    genes = set(gt['Gene1']).union(set(gt['Gene2']))
    n = len(genes)
    n_possible = n * (n - 1)

    if n_possible == 0:
        return float('nan'), float('nan'), float('nan')

    k_total = len(gt)
    has_type = 'Type' in gt.columns
    k_act = int((gt['Type'] == '+').sum()) if has_type else 0
    k_inh = int((gt['Type'] == '-').sum()) if has_type else 0

    b_total = k_total / n_possible if k_total > 0 else float('nan')
    b_act   = k_act   / n_possible if k_act   > 0 else float('nan')
    b_inh   = k_inh   / n_possible if k_inh   > 0 else float('nan')

    return b_total, b_act, b_inh


def _load_ratio_section(
    dataset_ids: List[str],
    dataset_paths: List[Path],
    csv_name: str,
    baselines: Dict[str, float],
) -> Dict[str, Dict[str, float]]:
    """
    Load median ratio values for one heatmap section.

    Reads {dataset_path}/{csv_name} for each dataset, divides each per-run
    value by the corresponding entry in baselines, then takes the median
    across runs. Datasets whose CSV is missing or whose baseline is nan or
    zero are skipped with a warning.

    Parameters
    ----------
    dataset_ids : list[str]
        Ordered dataset identifiers.
    dataset_paths : list[Path]
        Output directory for each dataset.
    csv_name : str
        Filename of the metric CSV (e.g. 'AUPRC.csv').
    baselines : dict[str, float]
        dataset_id -> random-predictor baseline for this metric.

    Returns
    -------
    dict[str, dict[str, float]]
        algo -> {dataset_id -> median ratio}, nan where unavailable.
    """
    result: Dict[str, Dict[str, float]] = {}

    for dataset_id, dataset_path in zip(dataset_ids, dataset_paths):
        csv_path = dataset_path / csv_name
        if not csv_path.exists():
            print(f"Warning: {csv_path} not found, skipping.")
            continue

        baseline = baselines.get(dataset_id, float('nan'))
        if math.isnan(baseline) or baseline == 0.0:
            print(
                f"Warning: baseline undefined for '{dataset_id}' "
                f"({csv_name}), skipping."
            )
            continue

        df = pd.read_csv(csv_path, index_col=0)
        for algo in df.index:
            vals = [v for v in df.loc[algo].tolist() if not math.isnan(v)]
            ratio = (
                statistics.median([v / baseline for v in vals])
                if vals else float('nan')
            )
            result.setdefault(str(algo), {})[dataset_id] = ratio

    return result


class PlotEPRHeatmap(Plotter):
    """
    Replicates Figure 4 of Pratapa et al. 2020 (BEELINE).

    Produces a four-section heatmap: median AUPRC ratio, median EPR ratio,
    median EPR for activating edges, and median EPR for inhibitory edges.
    Algorithms (rows) are sorted by decreasing median AUPRC ratio. Each
    section uses a distinct colour palette; cells below the random-predictor
    cutoff are shown in dark grey. Writes EPRSummary.pdf to the output
    directory.
    """

    def __call__(self, config: dict, output_dir: Path, root: Path) -> None:
        """
        Generate the EPR heatmap and write it to output_dir/EPRSummary.pdf.

        Parameters
        ----------
        config : dict
            Parsed YAML configuration.
        output_dir : Path
            Directory where EPRSummary.pdf is written.
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
            print("No datasets found for EPR heatmap.")
            return

        dataset_ids    = [r[0] for r in rows]
        # dataset_labels : list[str] — per-dataset display labels for column
        # headers; uses 'nickname' from config when set, else falls back to
        # dataset_id.
        dataset_labels = [r[1] for r in rows]
        dataset_paths  = [r[2] for r in rows]
        gt_paths       = [r[3] for r in rows]

        # Compute per-dataset random-predictor baselines for each metric type.
        b_total: Dict[str, float] = {}
        b_act:   Dict[str, float] = {}
        b_inh:   Dict[str, float] = {}

        for dataset_id, gt_path in zip(dataset_ids, gt_paths):
            if not gt_path.exists():
                print(f"Warning: {gt_path} not found, skipping {dataset_id}.")
                continue
            bt, ba, bi = _signed_random_baselines(gt_path)
            b_total[dataset_id] = bt
            b_act[dataset_id]   = ba
            b_inh[dataset_id]   = bi

        auprc_ratios    = _load_ratio_section(
            dataset_ids, dataset_paths, 'AUPRC.csv', b_total)
        # EarlyPrecision.csv already stores EPR values (ratio to random baseline),
        # so pass unit baselines to avoid dividing by the baseline a second time.
        unit_baselines  = {d: 1.0 for d in dataset_ids}
        epr_ratios      = _load_ratio_section(
            dataset_ids, dataset_paths, 'EarlyPrecision.csv', unit_baselines)
        epr_act_ratios  = _load_ratio_section(
            dataset_ids, dataset_paths, 'EarlyPrecisionActivation.csv', unit_baselines)
        epr_inh_ratios  = _load_ratio_section(
            dataset_ids, dataset_paths, 'EarlyPrecisionInhibitory.csv', unit_baselines)

        all_algos = sorted(
            set(auprc_ratios)
            | set(epr_ratios)
            | set(epr_act_ratios)
            | set(epr_inh_ratios)
        )
        if not all_algos:
            print("No data found for EPR heatmap.")
            return

        # Sort algorithms by decreasing median AUPRC ratio.
        def _median_auprc(algo: str) -> float:
            vals = [v for v in auprc_ratios.get(algo, {}).values()
                    if not math.isnan(v)]
            return statistics.median(vals) if vals else 0.0

        sorted_algos = sorted(all_algos, key=_median_auprc, reverse=True)
        n_algos    = len(sorted_algos)
        n_datasets = len(dataset_ids)

        # Build raw value arrays (n_algos × n_datasets).
        def _arr(section: Dict[str, Dict[str, float]]) -> np.ndarray:
            return np.array([
                [section.get(a, {}).get(d, float('nan')) for d in dataset_ids]
                for a in sorted_algos
            ])

        auprc_arr   = _arr(auprc_ratios)
        epr_arr     = _arr(epr_ratios)
        epr_act_arr = _arr(epr_act_ratios)
        epr_inh_arr = _arr(epr_inh_ratios)

        # --- Figure layout ---
        # 4 sections of n_datasets columns each, separated by gap columns.
        # Section x-starts: 1, n+2, 2n+3, 3n+4  (n = n_datasets)
        total_cols = n_datasets * 4 + 3
        pad        = 2
        height     = 7
        asp_ratio  = (total_cols + pad) / (n_algos + pad)
        fig_size   = (height * asp_ratio + 0.5, height)

        fig = plt.figure(figsize=(fig_size[0], fig_size[1] + 0.5))
        ax  = fig.add_subplot(111)

        _setup_heatmap_axes(ax, n_algos, sorted_algos, total_cols, pad)

        # Measure the maximum y-tick label width in axes-fraction coordinates so
        # the row backgrounds extend precisely to cover the algorithm labels.
        # draw() forces text layout before get_window_extent() is called.
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        max_label_px = max(
            (t.get_window_extent(renderer).width for t in ax.get_yticklabels()),
            default=0.0,
        )
        # Small padding (0.01) leaves a sliver of space between label and edge.
        label_frac = max_label_px / ax.get_window_extent(renderer).width + 0.01

        row_trans = blended_transform_factory(ax.transAxes, ax.transData)
        for row_idx in range(n_algos):
            bg = (0.9, 0.9, 0.9) if row_idx % 2 == 0 else (1.0, 1.0, 1.0)
            ax.add_artist(patches.Rectangle(
                (-label_frac, n_algos - row_idx - 0.5),
                width=1.0 + label_frac, height=1,
                transform=row_trans,
                clip_on=False,
                edgecolor=(1, 1, 1), facecolor=bg,
            ))

        # Color palettes — viridis for AUPRC, magma for the three EPR sections.
        # Magma is sampled from 20%–100% of its range so the low end is not
        # too dark; the full range (0%–100%) produces near-black at the bottom.
        auprc_palette = sns.color_palette("viridis", 11)
        epr_palette   = [mpl.cm.get_cmap("magma")(v) for v in np.linspace(0.2, 0.93, 11)]

        sections = [
            (auprc_arr,   auprc_palette, 'AUPRC Ratio',      1),
            (epr_arr,     epr_palette,   'EPR',              n_datasets + 2),
            (epr_act_arr, epr_palette,   'EPR (Activating)', n_datasets * 2 + 3),
            (epr_inh_arr, epr_palette,   'EPR (Inhibitory)', n_datasets * 3 + 4),
        ]

        for data, palette, label, col_x_start in sections:
            _draw_section(
                ax, data, n_algos, n_datasets,
                col_x_start=col_x_start,
                palette=palette,
                rand_cutoff=1.0,
                dataset_ids=dataset_labels,
                section_label=label,
            )

        # Legend geometry — positions are independent of column layout.
        # Three legends placed at the left third, center, and right third of
        # the axes width, each the same size.
        lw = 5.0 / (total_cols + 1)   # legend width as axes fraction
        lh = 0.5 / (n_algos + pad)    # legend height as axes fraction
        ly = -lh - 0.05               # just below the heatmap

        # Width that makes the random-predictor box square in figure inches.
        # lh is a fraction of axes height; multiply by height/width to get the
        # equivalent fraction of axes width.
        actual_fig_w = fig_size[0]
        actual_fig_h = fig_size[1] + 0.5
        rand_lw = lh * (actual_fig_h / actual_fig_w)

        for cx, palette, ticks, labels, w in [
            (1/6, auprc_palette, [0.5, len(auprc_palette) - 2], ['Low/Poor', 'High/Good'], lw),
            (1/2, None,          [0],                            ['Random Predictor'],      rand_lw),
            (5/6, epr_palette,   [0.5, len(epr_palette) - 2],   ['Low/Poor', 'High/Good'], lw),
        ]:
            legend_ax = ax.inset_axes([cx - w / 2, ly, w, lh])
            if palette is None:
                # Single grey cell for the random-predictor indicator.
                legend_ax.imshow(
                    np.zeros((1, 1)),
                    cmap=mpl.colors.ListedColormap([_DARK_COL]),
                    interpolation='nearest', aspect='auto',
                )
            else:
                legend_ax.imshow(
                    np.arange(len(palette)).reshape(1, len(palette)),
                    cmap=mpl.colors.ListedColormap(list(palette)),
                    interpolation='nearest', aspect='auto',
                )
            legend_ax.yaxis.set_ticks_position('none')
            legend_ax.xaxis.set_ticks_position('none')
            legend_ax.set_yticklabels([])
            legend_ax.set_xticks(ticks)
            legend_ax.set_xticklabels(labels, fontsize=12)

        stem = output_dir / 'EPRSummary'
        plt.savefig(stem.with_suffix('.pdf'), bbox_inches='tight')
        plt.savefig(stem.with_suffix('.png'), bbox_inches='tight')
        plt.close(fig)
        print(f"Saved EPR heatmap to {stem}.pdf and .png")
