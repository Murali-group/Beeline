"""
Shared drawing helpers for BEELINE heatmap plotters.

Provides the cell-drawing primitive (_flat_square), the section renderer
(_draw_section), and the axis-setup utility (_setup_heatmap_axes) used by
both PlotSummaryHeatmap and PlotEPRHeatmap.
"""
import math
from typing import List

import matplotlib.patches as patches
import numpy as np


# Color used for cells whose value falls below the random-predictor cutoff.
_DARK_COL = '#2E4053'


def _flat_square(ax, cx: float, cy: float, col) -> None:
    """
    Draw a fixed-size flat square cell centered at (cx, cy).

    Matches the 'hm' shape in old_plotter.py: a borderless Rectangle with a
    slight horizontal gap (width = 0.9) so adjacent cells do not merge visually.

    Parameters
    ----------
    ax : matplotlib Axes
    cx, cy : float
        Center coordinates in data units.
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
    show_all_values: bool = False,
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
    show_all_values : bool
        When True, annotate every cell (not just min/max) with its raw value.
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

        # Per-column min/max for normalization (all non-NaN values).
        col_vals = [data[r, col_idx] for r in range(n_algos)
                    if not math.isnan(data[r, col_idx])]
        col_min = min(col_vals) if col_vals else 0.0
        col_max = max(col_vals) if col_vals else 1.0

        # Minimum value that exceeds rand_cutoff — used for the low annotation
        # so the label always appears on a visible (non-dark) cell.
        visible_vals = [v for v in col_vals if v >= rand_cutoff]
        col_visible_min = min(visible_vals) if visible_vals else float('nan')

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

            # Column-wise min-max normalization for color selection.
            if col_max > col_min:
                norm_val = (raw - col_min) / (col_max - col_min)
            else:
                norm_val = 0.5

            col = palette[min(int(np.floor(norm_val * 10)), 10)]
            _flat_square(ax, cx, cy, col)

            # Annotate cells with the raw value.
            # Use one decimal place for values >= 1.0, two otherwise.
            fmt = f'{raw:.1f}' if raw >= 1.0 else f'{raw:.2f}'
            if show_all_values:
                text_col = text_col_max if norm_val >= 0.5 else text_col_min
                ax.text(cx, cy, fmt, fontsize=12,
                        ha='center', va='center', color=text_col,
                        bbox=dict(boxstyle='round', ec=(1,1,1,0), fc=(1,1,1,0)))
            elif norm_val >= 1.0:
                ax.text(cx, cy, fmt, fontsize=12,
                        ha='center', va='center', color=text_col_max,
                        bbox=dict(boxstyle='round', ec=(1,1,1,0), fc=(1,1,1,0)))
            elif not math.isnan(col_visible_min) and raw <= col_visible_min:
                ax.text(cx, cy, fmt, fontsize=12,
                        ha='center', va='center', color=text_col_min,
                        bbox=dict(boxstyle='round', ec=(1,1,1,0), fc=(1,1,1,0)))



def _setup_heatmap_axes(
    ax,
    n_algos: int,
    sorted_algos: List[str],
    total_cols: int,
    pad: int,
) -> None:
    """
    Configure tick labels, spines, and axis limits for a heatmap axes.

    Parameters
    ----------
    ax : matplotlib Axes
    n_algos : int
        Number of algorithm rows.
    sorted_algos : list[str]
        Algorithm names in draw order (best first); used for y-tick labels.
        Labels are reversed so the top row aligns with the best algorithm.
    total_cols : int
        Total number of data columns across all sections (used for xlim).
    pad : int
        Padding rows/columns added around the plot area (used for ylim).
    """
    ax.set_yticks(np.arange(0, n_algos + pad))
    # Labels run bottom-to-top; best algo (row 0, drawn at y=n_algos) must
    # align with the topmost label, so the list is reversed.
    ax.set_yticklabels([""] + sorted_algos[::-1] + [""], fontsize=12)

    ax.set_xticks(np.arange(0, total_cols + pad))
    ax.set_xticklabels([])
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_xlim(0, total_cols + 1)
    ax.set_ylim(0, n_algos + pad + 1)
