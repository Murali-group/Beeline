"""
PlotFigure5 — Reproduces the Figure 5 table layout from Pratapa et al. 2020.

Rows are (cell_type × network_type) dataset pairs derived by grouping the
500-gene and 1000-gene variants of each dataset. The left half of each row
shows data for TFs + 500 genes and the right half for TFs + 1000 genes.
Columns within each half: #TFs, #Genes, Density (plain text) followed by one
EPR column per algorithm from the config.

Datasets whose IDs do not contain a '-500-' or '-1000-' segment are skipped
with a warning. Algorithm column headers are truncated to four characters.
"""
import math
import re
import statistics
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.transforms import blended_transform_factory

from BLPlot._heatmap import _DARK_COL, _flat_square
from BLPlot.plotter import Plotter, get_algo_ids, iter_datasets_with_runs

plt.rcParams["font.size"] = 10

# Characters used for algorithm column header abbreviations.
_ALGO_ABBREV_LEN = 4

# Display labels for the three network-statistics columns.
_STAT_HEADERS: List[str] = ['#TFs', '#Genes', 'Density']
_N_STAT: int = len(_STAT_HEADERS)


def _row_key(dataset_id: str) -> str:
    """
    Return a grouping key for pairing 500- and 1000-gene dataset variants.

    Removes the gene-count segment from dataset_id so that e.g.
    'mHSC-E-500-string' and 'mHSC-E-1000-string' map to the same key
    ('mHSC-E-string').

    Parameters
    ----------
    dataset_id : str
        Dataset identifier (e.g. 'mHSC-E-500-string').

    Returns
    -------
    str
        Key with gene-count segment removed (e.g. 'mHSC-E-string').
    """
    if not isinstance(dataset_id, str):
        raise TypeError(f"dataset_id must be str, got {type(dataset_id)}")
    return re.sub(r'-(?:500|1000)-', '-', dataset_id)


def _gene_count(dataset_id: str) -> Optional[str]:
    """
    Return '500' or '1000' if the gene-count segment is present, else None.

    Parameters
    ----------
    dataset_id : str
        Dataset identifier.

    Returns
    -------
    str or None
        The gene-count string found in the identifier, or None.
    """
    if '-500-' in dataset_id:
        return '500'
    if '-1000-' in dataset_id:
        return '1000'
    return None


def _row_label(nickname: str) -> str:
    """
    Derive a row display label from a dataset nickname by removing the
    gene-count digit ('5' for 500, '1' for 1000).

    Nicknames follow the pattern {cell_abbr}{gene_digit}{net_abbr}
    (e.g. 'mHE5s'); removing the first digit yields the row label ('mHEs').
    Both 500- and 1000-gene variants of the same dataset produce the same
    label, so either may be passed.

    Parameters
    ----------
    nickname : str
        Dataset nickname as defined in the config (e.g. 'mHE5s').

    Returns
    -------
    str
        Nickname with the first digit character removed (e.g. 'mHEs').
    """
    if not isinstance(nickname, str):
        raise TypeError(f"nickname must be str, got {type(nickname)}")
    return re.sub(r'\d', '', nickname, count=1)


def _abbrev(algo_id: str, n: int = _ALGO_ABBREV_LEN) -> str:
    """
    Return the first n characters of algo_id as a column header abbreviation.

    Parameters
    ----------
    algo_id : str
        Algorithm identifier (e.g. 'GENIE3').
    n : int
        Maximum characters to keep (default: _ALGO_ABBREV_LEN).

    Returns
    -------
    str
        Truncated identifier (e.g. 'GENI').
    """
    if not isinstance(algo_id, str):
        raise TypeError(f"algo_id must be str, got {type(algo_id)}")
    return algo_id[:n]


def _compute_network_stats(gt_path: Path) -> Tuple[int, int, float]:
    """
    Compute summary statistics from a ground truth edge-list CSV.

    Self-loops (Gene1 == Gene2) are excluded. Gene1 is treated as the
    regulator (TF) column and Gene2 as the target column.

    Parameters
    ----------
    gt_path : Path
        Path to the ground truth CSV with columns Gene1 and Gene2.

    Returns
    -------
    tuple of (n_tfs, n_genes, density)
        n_tfs    — number of unique Gene1 (regulator) values.
        n_genes  — total unique genes across Gene1 ∪ Gene2.
        density  — k / (n_tfs * (n_genes - 1)) where k is the edge count;
                   nan when n_tfs * (n_genes - 1) == 0.
    """
    if not isinstance(gt_path, Path):
        raise TypeError(f"gt_path must be Path, got {type(gt_path)}")
    gt = pd.read_csv(gt_path, header=0)
    gt = gt[gt['Gene1'] != gt['Gene2']]
    tfs = set(gt['Gene1'].unique())
    genes = tfs | set(gt['Gene2'].unique())
    n_tfs = len(tfs)
    n_genes = len(genes)
    # Matches generateExpInputs.py: edges / (n_tfs * (n_genes - 1))
    denom = n_tfs * (n_genes - 1)
    density = len(gt) / denom if denom > 0 else float('nan')
    return n_tfs, n_genes, density


def _load_epr(dataset_path: Path, algo_ids: List[str]) -> Dict[str, float]:
    """
    Load median EPR values from EarlyPrecision.csv for the given algorithms.

    EarlyPrecision.csv has one row per algorithm (index column) and one data
    column per run. For single-run datasets there is exactly one value column.
    Only algorithms in algo_ids are returned; others are ignored.

    Parameters
    ----------
    dataset_path : Path
        Directory containing EarlyPrecision.csv.
    algo_ids : list[str]
        Algorithm identifiers to load.

    Returns
    -------
    dict[str, float]
        algo_id -> median EPR value. Absent when the CSV is missing or the
        algorithm is not present in the file.
    """
    if not isinstance(dataset_path, Path):
        raise TypeError(f"dataset_path must be Path, got {type(dataset_path)}")
    csv_path = dataset_path / 'EarlyPrecision.csv'
    if not csv_path.exists():
        return {}
    df = pd.read_csv(csv_path, index_col=0)
    result: Dict[str, float] = {}
    for algo in algo_ids:
        if algo not in df.index:
            continue
        vals = []
        for v in df.loc[algo].tolist():
            try:
                fv = float(v)
                if not math.isnan(fv):
                    vals.append(fv)
            except (TypeError, ValueError):
                pass
        if vals:
            result[algo] = statistics.median(vals)
    return result


def _draw_epr_col(
    ax,
    col_data: List[float],
    n_rows: int,
    cx: float,
    palette: list,
    rand_cutoff: float = 1.0,
    col_min: float = None,
    col_max: float = None,
    row_annotate: List[bool] = None,
) -> None:
    """
    Draw one EPR data column onto ax.

    Cells whose value falls below rand_cutoff are filled with _DARK_COL.
    All other cells are colored by linear normalization between col_min (or
    rand_cutoff when col_min is None) and col_max, mapped onto the palette.
    Value labels are drawn only on rows where row_annotate[row_idx] is True
    (or on all visible cells when row_annotate is None).

    Parameters
    ----------
    ax : matplotlib Axes
        Target axes.
    col_data : list[float]
        EPR value per row (index 0 = top row of the figure); nan = missing.
    n_rows : int
        Total number of data rows (sets y coordinate of row index 0).
    cx : float
        X centre of this column in data units.
    palette : list
        Color list; index 0 = low end, index -1 = high end.
    rand_cutoff : float
        Values strictly below this threshold are shown as _DARK_COL (default
        1.0, matching the EPR random-predictor baseline).
    col_min : float or None
        Lower bound for colour normalisation. When None, defaults to
        rand_cutoff. Pass a half-wide minimum so that all algorithm columns
        within the same half share the same colour scale floor.
    col_max : float or None
        Maximum value used for normalization. When None, defaults to the
        maximum non-NaN value in col_data. Pass a half-wide max so that all
        algorithm columns within the same half share the same colour scale.
    row_annotate : list[bool] or None
        Per-row flag controlling whether to draw the value label. When None,
        all visible (non-dark) cells are annotated.
    """
    col_vals = [v for v in col_data if not math.isnan(v)]
    if col_max is None:
        col_max = max(col_vals) if col_vals else 1.0
    norm_floor = col_min if col_min is not None else rand_cutoff
    n_pal = len(palette)

    for row_idx, raw in enumerate(col_data):
        cy = n_rows - row_idx

        if math.isnan(raw):
            ax.text(cx, cy, '-', fontsize=9, fontweight='medium',
                    ha='center', va='center', color=_DARK_COL)
            continue

        if raw < rand_cutoff:
            _flat_square(ax, cx, cy, _DARK_COL)
            continue

        # Normalize raw value from [norm_floor, col_max] to [0, 1].
        if col_max > norm_floor:
            norm_val = (raw - norm_floor) / (col_max - norm_floor)
        else:
            norm_val = 0.5

        col = palette[min(int(np.floor(norm_val * (n_pal - 1))), n_pal - 1)]
        _flat_square(ax, cx, cy, col)

        # Annotate only when this row is flagged (or always when no flags given).
        should_annotate = (row_annotate is None) or (row_idx < len(row_annotate) and row_annotate[row_idx])
        if should_annotate:
            fmt = f'{raw:.1f}' if raw >= 1.0 else f'{raw:.2f}'
            text_col = 'black' if norm_val >= 0.5 else 'white'
            ax.text(cx, cy, fmt, fontsize=9, fontweight='medium',
                    ha='center', va='center', color=text_col,
                    bbox=dict(boxstyle='round', ec=(1, 1, 1, 0), fc=(1, 1, 1, 0)))


class PlotFigure5(Plotter):
    """
    Reproduces the Figure 5 layout from Pratapa et al. 2020 (BEELINE).

    Datasets are paired by gene count: 500-gene variants occupy the left half
    of each row and 1000-gene variants the right half. Each half contains
    three network-statistics columns (#TFs, #Genes, Density) followed by one
    EPR column per algorithm. Missing data is shown as a dash or dark cell.
    Writes Figure5.pdf and Figure5.png to the output directory.
    """

    def __call__(self, config: dict, output_dir: Path, root: Path) -> None:
        """
        Generate the Figure 5 table and write it to output_dir/Figure5.pdf.

        Parameters
        ----------
        config : dict
            Parsed YAML configuration.
        output_dir : Path
            Directory where Figure5.pdf and Figure5.png are written.
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

        algo_ids = get_algo_ids(config)
        if not algo_ids:
            print("No algorithms found in config for Figure 5.")
            return

        # --- Build row map ---
        # row_map: OrderedDict[row_key, dict] — preserves config dataset order.
        # Each entry holds:
        #   'label'  : str — display label derived from the dataset nickname
        #   '500'    : (dataset_path, gt_path) or None
        #   '1000'   : (dataset_path, gt_path) or None
        row_map: 'OrderedDict[str, dict]' = OrderedDict()
        for dataset_id, nickname, dataset_path, gt_path, _ in iter_datasets_with_runs(config, root):
            gc = _gene_count(dataset_id)
            if gc is None:
                print(f"Warning: dataset '{dataset_id}' has no '-500-' or '-1000-' "
                      f"segment; skipping.")
                continue
            key = _row_key(dataset_id)
            if key not in row_map:
                row_map[key] = {'label': None, '500': None, '1000': None}
            row_map[key][gc] = (dataset_path, gt_path)
            # Use whichever variant we see first to set the label; both produce
            # the same output from _row_label since it removes the gene-count digit.
            if row_map[key]['label'] is None:
                row_map[key]['label'] = _row_label(nickname)

        rows = list(row_map.values())
        if not rows:
            print("No valid (500/1000-gene) datasets found for Figure 5.")
            return

        n_rows  = len(rows)
        n_algos = len(algo_ids)
        n_stat  = _N_STAT

        # --- Column x-coordinate layout (1-indexed data units) ---
        # Left half:   stat cols at 1..n_stat, algo cols at n_stat+1..n_stat+n_algos
        # Gap column:  n_stat + n_algos + 1
        # Right half:  stat cols at gap+1..gap+n_stat, algo cols at gap+n_stat+1..gap+n_stat+n_algos
        gap_x      = n_stat + n_algos              # right half begins at gap_x + 1
        total_cols = gap_x + n_stat + n_algos    # x of the last column

        # --- Header y-coordinate layout ---
        # Data rows occupy y = 1 (bottom) .. n_rows (top).
        # y = n_rows + 0.6 : column labels (#TFs, #Genes, Density, algo abbrevs)
        # y = n_rows + 2   : sub-headers ("Network Statistics" / "EPR")
        # y = n_rows + 3   : top-level half headers ("TFs + 500 genes" / "TFs + 1000 genes")
        y_col = n_rows + 0.6
        y_sub = n_rows + 1.7
        y_top = n_rows + 2.5
        pad   = 4   # y-units above data rows reserved for headers

        # --- Figure size ---
        # cell_size: figure inches per data unit; controls cell size and spacing.
        cell_size = 0.6
        fig_w = (total_cols + 3) * cell_size   # +3 for the row-label margin
        fig_h = (n_rows + pad + 1) * cell_size
        fig = plt.figure(figsize=(fig_w, fig_h))
        ax  = fig.add_subplot(111)

        # --- Axes limits and tick setup ---
        ax.set_xlim(0, total_cols + 1)
        ax.set_ylim(0, n_rows + pad + 1)

        # Y-ticks are hidden; labels are drawn manually as right-aligned text
        # just to the left of x=0.5 (the white vertical line) so they sit
        # flush against the left data boundary.
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.yaxis.set_ticks_position('none')

        ax.set_xticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        # --- Alternating row background stripes ---
        # Extend left past the y-axis just enough to cover the row-label text.
        # The stripe extension is kept narrow (0.08 axes fraction) to avoid
        # the stripe visibly jutting out beyond the label area.
        row_trans = blended_transform_factory(ax.transAxes, ax.transData)
        for row_idx in range(n_rows):
            bg = (0.9, 0.9, 0.9) if row_idx % 2 == 0 else (1.0, 1.0, 1.0)
            ax.add_artist(patches.Rectangle(
                (-0.08, n_rows - row_idx - 0.5),
                width=1.08, height=1.0,
                transform=row_trans,
                clip_on=False,
                edgecolor=(1, 1, 1), facecolor=bg,
            ))

        # Draw row labels right-aligned at x = 0.3 (just left of the white
        # vertical line at x = 0.5). Row index 0 maps to y = n_rows (top).
        for row_idx, row in enumerate(rows):
            cy = n_rows - row_idx
            ax.text(0.3, cy, row['label'] or '', fontsize=9, fontweight='medium',
                    ha='right', va='center', color='black')

        # --- EPR color palette: magma sampled from 20%–93% of its range ---
        # Matches the palette used in PlotEPRHeatmap.
        epr_palette = [mpl.cm.get_cmap("magma")(v)
                       for v in np.linspace(0.2, 0.93, 11)]

        # --- Phase 1: collect data for both halves before drawing ---
        # epr_by_half[gc][algo_id] = [float per row]
        # stat_by_half[gc] = [(n_tfs, n_genes, density) or None] per row
        # half_max_by_half[gc] / half_min_by_half[gc] = shared scale bounds per half
        epr_by_half:      Dict[str, Dict[str, List[float]]] = {}
        stat_by_half:     Dict[str, list]                   = {}
        half_max_by_half: Dict[str, float]                  = {}
        half_min_by_half: Dict[str, float]                  = {}

        for half_gc, x_offset in [('500', 0), ('1000', gap_x)]:
            epr_col_data: Dict[str, List[float]] = {a: [] for a in algo_ids}
            stat_col_data: list = []

            for row in rows:
                entry = row[half_gc]
                if entry is None:
                    stat_col_data.append(None)
                    for a in algo_ids:
                        epr_col_data[a].append(float('nan'))
                    continue

                dataset_path, gt_path = entry

                if gt_path.exists():
                    stat_col_data.append(_compute_network_stats(gt_path))
                else:
                    print(f"Warning: {gt_path} not found; skipping network stats.")
                    stat_col_data.append(None)

                epr_vals = _load_epr(dataset_path, algo_ids)
                for a in algo_ids:
                    epr_col_data[a].append(epr_vals.get(a, float('nan')))

            half_non_dark = [
                v for algo in algo_ids for v in epr_col_data[algo]
                if not math.isnan(v) and v >= 1.0
            ]
            epr_by_half[half_gc]      = epr_col_data
            stat_by_half[half_gc]     = stat_col_data
            half_max_by_half[half_gc] = max(half_non_dark) if half_non_dark else 1.0
            half_min_by_half[half_gc] = min(half_non_dark) if half_non_dark else 1.0

        # --- Phase 2 & 3: per-row max/min and annotation masks per half ---
        # Each half is treated independently: within a half, the highest and
        # lowest non-dark EPR values per row are annotated (up to 4 annotations
        # total across both halves when max ≠ min in each).
        row_annotate_map: Dict[str, Dict[str, List[bool]]] = {}
        for half_gc in ('500', '1000'):
            row_max = [float('-inf')] * n_rows
            row_min = [float('inf')]  * n_rows
            for algo in algo_ids:
                for r, v in enumerate(epr_by_half[half_gc][algo]):
                    if not math.isnan(v) and v >= 1.0:
                        if v > row_max[r]:
                            row_max[r] = v
                        if v < row_min[r]:
                            row_min[r] = v

            row_annotate_map[half_gc] = {}
            for algo in algo_ids:
                flags: List[bool] = []
                for r, v in enumerate(epr_by_half[half_gc][algo]):
                    if math.isnan(v) or v < 1.0:
                        flags.append(False)
                    else:
                        is_max = math.isclose(v, row_max[r], rel_tol=1e-9)
                        is_min = math.isclose(v, row_min[r], rel_tol=1e-9)
                        flags.append(is_max or is_min)
                row_annotate_map[half_gc][algo] = flags

        # --- Phase 4: draw both halves ---
        for half_gc, x_offset in [('500', 0), ('1000', gap_x)]:
            stat_x = [x_offset + i + 1 for i in range(n_stat)]
            algo_x = [x_offset + n_stat + i + 1 for i in range(n_algos)]

            epr_col_data = epr_by_half[half_gc]
            stat_col_data = stat_by_half[half_gc]
            half_max = half_max_by_half[half_gc]
            half_min = half_min_by_half[half_gc]

            # Draw network-statistics text columns.
            for row_idx, stats in enumerate(stat_col_data):
                if stats is None:
                    continue
                n_tfs, n_genes, density = stats
                cy = n_rows - row_idx
                for cx, text in zip(stat_x, [str(n_tfs), str(n_genes), f'{density:.2f}']):
                    ax.text(cx, cy, text, fontsize=9, fontweight='medium',
                            ha='center', va='center', color='black')

            # Draw EPR columns with per-row annotation masks.
            for algo_id, cx in zip(algo_ids, algo_x):
                _draw_epr_col(
                    ax, epr_col_data[algo_id], n_rows, cx,
                    epr_palette, rand_cutoff=1.0, col_min=half_min, col_max=half_max,
                    row_annotate=row_annotate_map[half_gc][algo_id],
                )

            # --- Column-level headers (#TFs, #Genes, Density, algo abbrevs) ---
            for cx, hdr in zip(stat_x, _STAT_HEADERS):
                ax.text(cx - 0.3, y_col, hdr, fontsize=8, fontweight='medium',
                        ha='left', va='bottom', rotation=30,
                        bbox=dict(boxstyle='round', ec=(1, 1, 1, 0), fc=(1, 1, 1, 0)))
            for algo_id, cx in zip(algo_ids, algo_x):
                ax.text(cx - 0.3, y_col, _abbrev(algo_id), fontsize=8, fontweight='medium',
                        ha='left', va='bottom', rotation=30,
                        bbox=dict(boxstyle='round', ec=(1, 1, 1, 0), fc=(1, 1, 1, 0)))

            # --- Sub-headers: "Network Statistics" and "EPR" ---
            stat_center = x_offset + (n_stat + 1) / 2
            epr_center  = x_offset + n_stat + (n_algos + 1) / 2
            ax.text(stat_center, y_sub, 'Network Statistics', fontsize=9, fontweight='semibold',
                    ha='center', va='center',
                    bbox=dict(boxstyle='round', ec=(1, 1, 1), fc=(1, 1, 1)))
            ax.text(epr_center, y_sub, 'EPR', fontsize=9, fontweight='semibold',
                    ha='center', va='center',
                    bbox=dict(boxstyle='round', ec=(1, 1, 1), fc=(1, 1, 1)))

            # Separator line between top-level header and sub-headers.
            # The two halves share the boundary x = n_stat+n_algos+0.5; a small
            # inset (0.3 data units) on each side leaves a visible gap between
            # the left-half line end and the right-half line start.
            is_left = (x_offset == 0)
            x0_line = x_offset + 0.5 + (0.0 if is_left else 0.3)
            x1_line = x_offset + n_stat + n_algos + 0.5 - (0.3 if is_left else 0.0)
            ax.plot([x0_line, x1_line], [y_sub + 0.5, y_sub + 0.5],
                    color='black', linewidth=2.0, solid_capstyle='butt')

            # --- Top-level half header: "TFs + N genes" ---
            half_center = x_offset + (n_stat + n_algos) / 2
            ax.text(half_center, y_top, f'TFs + {half_gc} genes',
                    fontsize=10, ha='center', va='center', fontweight='bold',
                    bbox=dict(boxstyle='round', ec=(1, 1, 1), fc=(1, 1, 1)))


        # --- Thin white line at the left edge of the 500-gene data area ---
        # Masks the sliver of alternating row stripe that bleeds past the data boundary.
        ax.axvline(x=0.5, color='white', linewidth=1.5, zorder=3)

        # --- Legend ---
        lw = 5.0 / (total_cols + 1)   # legend width as fraction of axes width
        lh = 0.5 / (n_rows + pad)     # legend height as fraction of axes height
        ly = -lh - 0.06               # y position: just below the heatmap

        # Width that makes the random-predictor box square in figure inches.
        rand_lw = lh * (fig_h / fig_w)

        # Centre of the left and right halves in axes fraction.
        left_half_cx  = (gap_x / 2) / (total_cols + 1)
        right_half_cx = (gap_x + (total_cols - gap_x) / 2) / (total_cols + 1)

        for cx, palette_arg, ticks, tick_labels, w in [
            (left_half_cx,  None,       [0],                          ['Random Predictor'],      rand_lw),
            (right_half_cx, epr_palette, [0.5, len(epr_palette) - 2], ['Low/Poor', 'High/Good'], lw),
        ]:
            legend_ax = ax.inset_axes([cx - w / 2, ly, w, lh])
            if palette_arg is None:
                legend_ax.imshow(
                    np.zeros((1, 1)),
                    cmap=mpl.colors.ListedColormap([_DARK_COL]),
                    interpolation='nearest', aspect='auto',
                )
            else:
                legend_ax.imshow(
                    np.arange(len(palette_arg)).reshape(1, len(palette_arg)),
                    cmap=mpl.colors.ListedColormap(list(palette_arg)),
                    interpolation='nearest', aspect='auto',
                )
            legend_ax.yaxis.set_ticks_position('none')
            legend_ax.xaxis.set_ticks_position('none')
            legend_ax.set_yticklabels([])
            legend_ax.set_xticks(ticks)
            legend_ax.set_xticklabels(tick_labels, fontsize=9, fontweight='semibold')

        output_dir.mkdir(parents=True, exist_ok=True)
        stem = output_dir / 'Figure5'
        plt.savefig(stem.with_suffix('.pdf'), bbox_inches='tight')
        plt.savefig(stem.with_suffix('.png'), bbox_inches='tight')
        plt.close(fig)
        print(f"Saved Figure 5 to {stem}.pdf and .png")
