"""
PlotFigure6 — Reproduces the Figure 6 table layout from Pratapa et al. 2020.

All data is supplied via a standalone YAML config file; no algorithm output
files are read. Each row is one algorithm. Columns are split into four sections:
  Properties  — plain text (Category, Addl. Inputs) and boolean icon columns
  Accuracy    — heatmap with a viridis palette (Synthetic, Curated, scRNA-Seq)
  Stability   — heatmap with a purple palette (Datasets, Runs, Dropouts, Pseudotime)
  Scalability — heatmap + text labels, purple palette, split into Time and Memory
                sub-sections each with 100 / 500 / 1k / 2k gene columns

Heatmap cells are normalised per column across all non-missing algorithms.
Scalability values are auto-parsed from strings (e.g. "1m", "0.5G"); lower
time / memory maps to a higher (lighter / better) colour.

Missing data can be represented either by omitting the key or by setting it to
null in the YAML; both render as a dash with no coloured background.

Writes Figure6.pdf and Figure6.png to the output directory.
"""
import logging
import math
import textwrap
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

plt.rcParams["font.size"] = 10
plt.rcParams["font.family"] = "Nimbus Sans"

# Suppress "can not be subsetted into a Type 3 font" log messages from the
# PDF backend when embedding OTF fonts that lack a subsetting table.
logging.getLogger('matplotlib.backends.backend_pdf').setLevel(logging.ERROR)

# ─── Symbols ──────────────────────────────────────────────────────────────────
_CHECK     = '✓'
_CROSS     = '✗'
_DASH      = '–'
_CHECK_COL = '#27AE60'   # green for true
_CROSS_COL = '#E74C3C'   # red for false

# ─── Column definitions ───────────────────────────────────────────────────────
# Each entry: (config_key, display_label, width_in_data_units)

# Properties section: text and boolean columns
_PROP_COLS: List[Tuple[str, str, float]] = [
    ('category',          'Category',      1.8),
    ('additional_inputs', 'Addl. Inputs',  1.6),
    ('time_ordered',      'Time ordered?', 1.0),
    ('directed',          'Directed?',     1.0),
    ('signed',            'Signed?',       1.0),
]

# Accuracy section: heatmap columns (viridis palette, higher = better)
_ACC_COLS: List[Tuple[str, str, float]] = [
    ('synthetic', 'Synthetic', 1.0),
    ('curated',   'Curated',   1.0),
    ('scrna_seq', 'scRNA-Seq', 1.0),
]

# Stability section: heatmap columns (purple palette, higher = better)
_STAB_COLS: List[Tuple[str, str, float]] = [
    ('datasets',   'Datasets',   1.0),
    ('runs',       'Runs',       1.0),
    ('dropouts',   'Dropouts',   1.0),
    ('pseudotime', 'Pseudotime', 1.0),
]

# Scalability sub-section gene-count keys (shared by Time and Memory)
_SCALE_KEYS: List[str] = ['100', '500', '1k', '2k']
# Width per scalability column — slightly wider than heatmap-only cols to
# accommodate text labels such as "30m" or "0.5G".
_SCALE_COL_W: float = 1.3

# Gap (data units) between sections / sub-sections
_SEC_GAP:     float = 0.5   # between Properties, Accuracy, Stability, Scalability
_SUBSEC_GAP:  float = 0.3   # between Time and Memory within Scalability

# Palette size (number of discrete colour steps)
_N_PAL: int = 11


# ─── Palette helpers ──────────────────────────────────────────────────────────

def _acc_palette() -> list:
    """
    Viridis-sampled palette for Accuracy columns.
    Low score → dark blue-green; high score → bright yellow-green.

    Returns
    -------
    list of RGBA tuples, length _N_PAL.
    """
    return [mpl.cm.get_cmap("viridis")(v) for v in np.linspace(0.1, 0.95, _N_PAL)]


def _stab_palette() -> list:
    """
    Cubehelix palette (reversed) for Stability columns, matching the
    Summary heatmap (sns.cubehelix_palette(11, reverse=True)).
    Low score → dark; high score → light.

    Returns
    -------
    list of RGBA tuples, length _N_PAL.
    """
    return list(sns.cubehelix_palette(_N_PAL, reverse=True))


def _scale_palette() -> list:
    """
    Magma-sampled palette for Scalability columns (Time and Memory).
    Low normalised score (slow/large) → dark; high (fast/small) → light.
    Sampled from 20%–93% of the magma range, matching the EPR summary plot.

    Returns
    -------
    list of RGBA tuples, length _N_PAL.
    """
    return [mpl.cm.get_cmap("magma")(v) for v in np.linspace(0.2, 0.93, _N_PAL)]


# ─── Column layout ────────────────────────────────────────────────────────────

def _col_layout(
    col_defs: List[Tuple[str, str, float]],
    x_start: float,
) -> Tuple[List[Tuple[str, str, float, float]], float]:
    """
    Compute centre x-coordinates for a list of column definitions.

    Parameters
    ----------
    col_defs : list of (key, label, width)
        Column definitions in left-to-right order.
    x_start : float
        Left edge of the first column in data units.

    Returns
    -------
    layout : list of (key, label, width, center_x)
    x_end  : float — right edge of the last column.
    """
    layout = []
    x = x_start
    for key, label, w in col_defs:
        layout.append((key, label, w, x + w / 2))
        x += w
    return layout, x


# ─── Parsing helpers ──────────────────────────────────────────────────────────

def _parse_time_s(s: Optional[str]) -> float:
    """
    Parse a time label string to seconds.

    Recognised suffixes: s (seconds), m (minutes), h (hours), d (days).
    A leading '>' is accepted and treated as 1.5× the stated value to
    position '>X' cells as worse-than but finite.
    Returns nan for None, missing, or unrecognised strings.

    Parameters
    ----------
    s : str or None
        Time label such as '1s', '30m', '1h', '>1d'.

    Returns
    -------
    float
        Parsed duration in seconds, or nan.
    """
    if s is None:
        return float('nan')
    s = str(s).strip()
    if s in ('', _DASH, '-'):
        return float('nan')
    mult = 1.5 if s.startswith('>') else 1.0
    s = s.lstrip('>')
    try:
        if s.endswith('s'):
            return float(s[:-1]) * mult
        if s.endswith('m'):
            return float(s[:-1]) * 60 * mult
        if s.endswith('h'):
            return float(s[:-1]) * 3600 * mult
        if s.endswith('d'):
            return float(s[:-1]) * 86400 * mult
    except ValueError:
        pass
    return float('nan')


def _parse_memory_gb(s: Optional[str]) -> float:
    """
    Parse a memory label string to gigabytes.

    Recognised suffixes: M (megabytes), G (gigabytes).
    A leading '>' is accepted and treated as 1.5× the stated value.
    Returns nan for None, missing, or unrecognised strings.

    Parameters
    ----------
    s : str or None
        Memory label such as '0.1G', '1G', '1M', '>4G'.

    Returns
    -------
    float
        Parsed size in GB, or nan.
    """
    if s is None:
        return float('nan')
    s = str(s).strip()
    if s in ('', _DASH, '-'):
        return float('nan')
    mult = 1.5 if s.startswith('>') else 1.0
    s = s.lstrip('>')
    try:
        if s.endswith('G'):
            return float(s[:-1]) * mult
        if s.endswith('M'):
            return float(s[:-1]) / 1024 * mult
    except ValueError:
        pass
    return float('nan')


# ─── Normalisation helpers ────────────────────────────────────────────────────

def _normalize(values: List[float]) -> List[float]:
    """
    Normalise a list of values to [0, 1]; nan values remain nan.

    When all non-nan values are equal, each is mapped to 0.5.

    Parameters
    ----------
    values : list of float
        Raw scores; nan indicates a missing cell.

    Returns
    -------
    list of float
        Normalised scores in [0, 1] or nan.
    """
    valid = [v for v in values if not math.isnan(v)]
    if not valid:
        return [float('nan')] * len(values)
    lo, hi = min(valid), max(valid)
    if lo == hi:
        return [float('nan') if math.isnan(v) else 0.5 for v in values]
    return [float('nan') if math.isnan(v) else (v - lo) / (hi - lo) for v in values]


def _normalize_rank_inverted(values: List[float]) -> List[float]:
    """
    Rank-based normalisation, inverted so the smallest value maps to 1.0
    (best / lightest) and the largest maps to 0.0 (worst / darkest).

    Each distinct non-nan value is assigned an equal step on the colour
    scale regardless of its numeric distance from its neighbours — e.g.
    [1, 2, 5, 1000000] all occupy evenly spaced steps. Ties share the
    same colour. nan values remain nan.

    Parameters
    ----------
    values : list of float
        Raw scalability values (seconds or GB) across all columns and
        algorithms (flattened).

    Returns
    -------
    list of float
        Rank-based normalised scores in [0, 1] or nan.
    """
    unique_sorted = sorted(set(v for v in values if not math.isnan(v)))
    n = len(unique_sorted)
    if n == 0:
        return [float('nan')] * len(values)
    if n == 1:
        return [float('nan') if math.isnan(v) else 0.5 for v in values]
    rank_map = {v: i / (n - 1) for i, v in enumerate(unique_sorted)}
    # Invert: smallest rank → 1.0 (best), largest → 0.0 (worst)
    return [float('nan') if math.isnan(v) else 1.0 - rank_map[v] for v in values]


# ─── Drawing helpers ──────────────────────────────────────────────────────────

def _flat_rect(ax, cx: float, cy: float, w: float, color) -> None:
    """
    Draw a filled rectangle centred at (cx, cy) with a small horizontal gap
    so adjacent cells do not merge visually. Width is rendered at 90% of the
    column slot width (matching the _flat_square convention in _heatmap.py).

    Parameters
    ----------
    ax : matplotlib Axes
    cx : float — x centre in data units.
    cy : float — y centre in data units.
    w  : float — column slot width in data units (rendered width = 0.9 * w).
    color      — fill colour (any matplotlib colour spec).
    """
    draw_w = w * 0.93
    ax.add_patch(patches.Rectangle(
        (cx - draw_w / 2, cy - 0.5), draw_w, 1.0,
        linewidth=0, facecolor=color, zorder=2,
    ))


def _palette_color(norm_val: float, palette: list):
    """
    Map a normalised value in [0, 1] to a colour from the palette list.

    Parameters
    ----------
    norm_val : float — normalised score in [0, 1].
    palette  : list  — list of RGBA tuples, length _N_PAL.

    Returns
    -------
    RGBA tuple.
    """
    idx = min(int(np.floor(norm_val * (len(palette) - 1))), len(palette) - 1)
    return palette[idx]


# ─── Main class ───────────────────────────────────────────────────────────────

class PlotFigure6:
    """
    Generates the Figure 6 summary table from a standalone YAML config.

    All data (algorithm properties, accuracy scores, stability scores,
    scalability time/memory labels) is read from the config dict. No
    algorithm output files are accessed.

    The config must have a top-level 'algorithms' list; each entry is a
    dict with the keys described in the module docstring.
    """

    def __call__(self, config: dict, output_dir: Path) -> None:
        """
        Generate Figure 6 and write it to output_dir/Figure6.pdf/.png.

        Parameters
        ----------
        config : dict
            Parsed Figure 6 YAML configuration.
        output_dir : Path
            Directory where Figure6.pdf and Figure6.png are written.

        Returns
        -------
        None
        """
        if not isinstance(config, dict):
            raise TypeError(f"config must be dict, got {type(config)}")
        if not isinstance(output_dir, Path):
            raise TypeError(f"output_dir must be Path, got {type(output_dir)}")

        algorithms = config.get('algorithms', [])
        if not algorithms:
            print("No algorithms found in Figure 6 config.")
            return

        n_rows    = len(algorithms)
        acc_pal   = _acc_palette()
        stab_pal  = _stab_palette()
        scale_pal = _scale_palette()

        # ── Build column layouts ──────────────────────────────────────────────
        x = 0.5   # left edge of first column
        prop_layout,  x = _col_layout(_PROP_COLS, x)
        x += _SEC_GAP
        acc_layout,   x = _col_layout(_ACC_COLS, x)
        x += _SEC_GAP
        stab_layout,  x = _col_layout(_STAB_COLS, x)
        x += _SEC_GAP
        scale_col_defs = [(k, k, _SCALE_COL_W) for k in _SCALE_KEYS]
        time_layout,  x = _col_layout(scale_col_defs, x)
        x += _SUBSEC_GAP
        mem_layout,   x = _col_layout(scale_col_defs, x)
        total_width = x + 0.5   # right margin

        # ── Pre-compute normalised scores per column ──────────────────────────
        def _raw(algo: dict, section: str, key: str) -> float:
            # Return float raw score or nan for missing / null.
            sec = algo.get(section) or {}
            v = sec.get(key)
            if v is None:
                return float('nan')
            try:
                f = float(v)
                return float('nan') if math.isnan(f) else f
            except (TypeError, ValueError):
                return float('nan')

        def _scale_label(algo: dict, sub: str, key: str) -> Optional[str]:
            # Return the string label for a scalability cell or None.
            scale = algo.get('scalability') or {}
            sub_sec = scale.get(sub) or {}
            v = sub_sec.get(key)
            return str(v) if v is not None else None

        # Accuracy: per-column normalisation, higher = better
        acc_norm: Dict[str, List[float]] = {}
        for key, _, _ in _ACC_COLS:
            raw = [_raw(a, 'accuracy', key) for a in algorithms]
            acc_norm[key] = _normalize(raw)

        # Stability: per-column normalisation, higher = better
        stab_norm: Dict[str, List[float]] = {}
        for key, _, _ in _STAB_COLS:
            raw = [_raw(a, 'stability', key) for a in algorithms]
            stab_norm[key] = _normalize(raw)

        # Scalability: parse labels, collect all raw values per metric, then
        # normalise globally across all gene-count columns so the colour scale
        # is consistent (e.g. "1m for 100 genes" and "1m for 2k genes" get the
        # same colour).  Invert so lower time / memory = higher (lighter) colour.
        time_norm:   Dict[str, List[float]]          = {}
        time_labels: Dict[str, List[Optional[str]]]  = {}
        mem_norm:    Dict[str, List[float]]           = {}
        mem_labels:  Dict[str, List[Optional[str]]]  = {}

        # Collect raw values per key (preserving order for later reshaping).
        _time_raw: Dict[str, List[float]] = {}
        _mem_raw:  Dict[str, List[float]] = {}
        for k in _SCALE_KEYS:
            t_lbls = [_scale_label(a, 'time',   k) for a in algorithms]
            m_lbls = [_scale_label(a, 'memory', k) for a in algorithms]
            time_labels[k] = t_lbls
            mem_labels[k]  = m_lbls
            _time_raw[k] = [_parse_time_s(s)    for s in t_lbls]
            _mem_raw[k]  = [_parse_memory_gb(s) for s in m_lbls]

        # Flatten all raw values and normalise with fixed log2 bounds so the
        # colour position of a given time/memory value is independent of what
        # other algorithms are present in the config.
        n_algo = len(algorithms)
        flat_time_norm = _normalize_rank_inverted(
            [v for k in _SCALE_KEYS for v in _time_raw[k]]
        )
        flat_mem_norm = _normalize_rank_inverted(
            [v for k in _SCALE_KEYS for v in _mem_raw[k]]
        )
        for i, k in enumerate(_SCALE_KEYS):
            time_norm[k] = flat_time_norm[i * n_algo : (i + 1) * n_algo]
            mem_norm[k]  = flat_mem_norm[ i * n_algo : (i + 1) * n_algo]

        # ── Figure and axes ───────────────────────────────────────────────────
        pad       = 4      # y-units above data reserved for headers
        cell_size = 0.5    # figure inches per data unit

        # Extra left margin for row labels
        label_margin = 1.5
        fig_w = (total_width + label_margin) * cell_size
        leg_rows = 2    # extra rows below data for legends
        fig_h = (n_rows + pad + 1 + leg_rows) * cell_size

        fig = plt.figure(figsize=(fig_w, fig_h))
        ax  = fig.add_subplot(111)

        ax.set_xlim(0, total_width)
        ax.set_ylim(-leg_rows, n_rows + pad + 1)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        # ── Header y-positions ────────────────────────────────────────────────
        y_col = n_rows + 0.7    # column labels (rotated)
        y_sub = n_rows + 1.8    # sub-section headers (Time / Memory)
        y_sec = n_rows + 2.8    # section headers (Properties, Accuracy, …)

        # ── Alternating row backgrounds ───────────────────────────────────────
        # Extends left to cover row label area.
        row_trans = mpl.transforms.blended_transform_factory(
            ax.transAxes, ax.transData,
        )
        for row_idx in range(n_rows):
            bg = (0.9, 0.9, 0.9) if row_idx % 2 == 0 else (1.0, 1.0, 1.0)
            ax.add_artist(patches.Rectangle(
                (-0.08, n_rows - row_idx - 0.5),
                width=1.08, height=1.0,
                transform=row_trans, clip_on=False,
                edgecolor=(1, 1, 1), facecolor=bg,
            ))

        # ── Row labels (algorithm names) ──────────────────────────────────────
        label_x = 0.5 - 0.15   # just left of data columns
        for row_idx, algo in enumerate(algorithms):
            cy = n_rows - row_idx
            ax.text(label_x, cy, algo.get('algorithm_id', ''),
                    fontsize=9, fontweight='medium',
                    ha='right', va='center', color='black')

        # ── Draw data cells ───────────────────────────────────────────────────
        for row_idx, algo in enumerate(algorithms):
            cy = n_rows - row_idx

            # Properties: text columns
            for key, label, w, cx in prop_layout:
                if key == 'category':
                    raw = algo.get('category')
                    val = _DASH if (raw is None or str(raw).strip() in ('None', 'null', 'NULL', '')) else str(raw)
                    is_dash = val == _DASH
                    ax.text(cx if is_dash else cx - w * 0.45,
                            cy, val, fontsize=8,
                            fontweight='bold' if is_dash else 'medium',
                            ha='center' if is_dash else 'left',
                            va='center', color='black')

                elif key == 'additional_inputs':
                    raw = algo.get('additional_inputs')
                    val = _DASH if (raw is None or str(raw).strip() in ('None', 'null', 'NULL', '')) else str(raw)
                    is_dash = val == _DASH
                    # Wrap long strings to at most 2 lines at ~12 chars per line.
                    display = val if is_dash else textwrap.fill(val, width=12, max_lines=3)
                    ax.text(cx if is_dash else cx - w * 0.45,
                            cy, display, fontsize=7,
                            fontweight='bold' if is_dash else 'medium',
                            ha='center' if is_dash else 'left',
                            va='center', color='black',
                            multialignment='left', linespacing=1.2)

                else:  # boolean column
                    val = algo.get(key)
                    if val is None:
                        ax.text(cx, cy, _DASH, fontsize=9, fontweight='bold',
                                ha='center', va='center', color='black')
                    elif bool(val):
                        ax.text(cx, cy, _CHECK, fontsize=10, fontweight='medium',
                                ha='center', va='center', color=_CHECK_COL,
                                fontfamily='DejaVu Sans')
                    else:
                        ax.text(cx, cy, _CROSS, fontsize=10, fontweight='medium',
                                ha='center', va='center', color=_CROSS_COL,
                                fontfamily='DejaVu Sans')

            # Accuracy: heatmap
            for key, label, w, cx in acc_layout:
                nv = acc_norm[key][row_idx]
                if math.isnan(nv):
                    ax.text(cx, cy, _DASH, fontsize=9, fontweight='bold',
                            ha='center', va='center', color='grey')
                else:
                    _flat_rect(ax, cx, cy, w, _palette_color(nv, acc_pal))

            # Stability: heatmap
            for key, label, w, cx in stab_layout:
                nv = stab_norm[key][row_idx]
                if math.isnan(nv):
                    ax.text(cx, cy, _DASH, fontsize=9, fontweight='bold',
                            ha='center', va='center', color='grey')
                else:
                    _flat_rect(ax, cx, cy, w, _palette_color(nv, stab_pal))

            # Scalability Time: heatmap + text label
            for key, _, w, cx in time_layout:
                nv  = time_norm[key][row_idx]
                lbl = time_labels[key][row_idx]
                if lbl is None or math.isnan(nv):
                    ax.text(cx, cy, _DASH, fontsize=8, fontweight='bold',
                            ha='center', va='center', color='grey')
                else:
                    _flat_rect(ax, cx, cy, w, _palette_color(nv, scale_pal))
                    text_col = 'black' if nv >= 0.5 else 'white'
                    ax.text(cx, cy, lbl, fontsize=7, fontweight='medium',
                            ha='center', va='center', color=text_col)

            # Scalability Memory: heatmap + text label
            for key, _, w, cx in mem_layout:
                nv  = mem_norm[key][row_idx]
                lbl = mem_labels[key][row_idx]
                if lbl is None or math.isnan(nv):
                    ax.text(cx, cy, _DASH, fontsize=8, fontweight='bold',
                            ha='center', va='center', color='grey')
                else:
                    _flat_rect(ax, cx, cy, w, _palette_color(nv, scale_pal))
                    text_col = 'black' if nv >= 0.5 else 'white'
                    ax.text(cx, cy, lbl, fontsize=7, fontweight='medium',
                            ha='center', va='center', color=text_col)

        # ── Column labels (rotated) ───────────────────────────────────────────
        all_col_layouts = (
            prop_layout + acc_layout + stab_layout + time_layout + mem_layout
        )
        for key, label, w, cx in all_col_layouts:
            ax.text(cx - 0.15, y_col, label,
                    fontsize=8, fontweight='medium',
                    ha='left', va='bottom', rotation=45,
                    bbox=dict(boxstyle='round', ec=(1,1,1,0), fc=(1,1,1,0)))

        # ── Sub-section headers (Time / Memory) ──────────────────────────────
        time_cx = time_layout[0][3] + (time_layout[-1][3] - time_layout[0][3]) / 2
        mem_cx  = mem_layout[0][3]  + (mem_layout[-1][3]  - mem_layout[0][3])  / 2
        for cx, lbl in [(time_cx, 'Time'), (mem_cx, 'Memory')]:
            ax.text(cx, y_sub, lbl, fontsize=9, fontweight='medium',
                    ha='center', va='center',
                    bbox=dict(boxstyle='round', ec=(1,1,1), fc=(1,1,1)))

        # ── Section headers ───────────────────────────────────────────────────
        def _section_cx(layout):
            # Centre x of a section given its column layout list.
            return (layout[0][3] + layout[-1][3]) / 2

        prop_cx  = _section_cx(prop_layout)
        acc_cx   = _section_cx(acc_layout)
        stab_cx  = _section_cx(stab_layout)
        scale_cx = (time_layout[0][3] + mem_layout[-1][3]) / 2

        for cx, lbl in [
            (prop_cx,  'Properties'),
            (acc_cx,   'Accuracy'),
            (stab_cx,  'Stability'),
            (scale_cx, 'Scalability (genes)'),
        ]:
            ax.text(cx, y_sec, lbl, fontsize=10, fontweight='bold',
                    ha='center', va='center',
                    bbox=dict(boxstyle='round', ec=(1,1,1), fc=(1,1,1)))

        # ── Thin white line at the left edge of the data area ────────────────
        # Masks any row-background bleed past the left data boundary.
        prop_x0 = prop_layout[0][3] - prop_layout[0][2] / 2
        ax.axvline(x=prop_x0, color='white', linewidth=1.5, zorder=3)

        # ── Separator lines ───────────────────────────────────────────────────
        # Draw one segment per section with gaps at the section boundaries.
        def _sec_x(layout_left, layout_right):
            # (x_left_edge, x_right_edge) for a section spanning two layouts.
            x0 = layout_left[0][3]  - layout_left[0][2]  / 2
            x1 = layout_right[-1][3] + layout_right[-1][2] / 2
            return x0, x1

        sections_x = [
            _sec_x(prop_layout,  prop_layout),
            _sec_x(acc_layout,   acc_layout),
            _sec_x(stab_layout,  stab_layout),
            _sec_x(time_layout,  mem_layout),   # Scalability spans both sub-sections
        ]

        # Line under section headers only (y_sub + 0.5); the top-of-data line
        # (n_rows + 0.5) is removed per-request.
        for x0, x1 in sections_x:
            ax.plot(
                [x0, x1], [y_sub + 0.5, y_sub + 0.5],
                color='black', linewidth=1.5, solid_capstyle='butt',
            )

        # Thin lines under "Time" and "Memory" sub-headers only.
        time_x0, time_x1 = _sec_x(time_layout, time_layout)
        mem_x0,  mem_x1  = _sec_x(mem_layout,  mem_layout)
        inset = 0.25   # horizontal inset from section edges
        for x0, x1 in [(time_x0, time_x1), (mem_x0, mem_x1)]:
            ax.plot(
                [x0 + inset, x1 - inset], [y_sub - 0.35, y_sub - 0.35],
                color='black', linewidth=0.8, solid_capstyle='butt',
            )

        # ── Legends ───────────────────────────────────────────────────────────
        # Matches PlotSummaryHeatmap: inset_axes + imshow strip, tick labels for
        # Low/Poor and High/Good. Spines hidden to avoid black outline.
        y_span = n_rows + leg_rows + pad + 1
        lw = 5.0 / total_width          # legend width as axes fraction
        lh = 0.5 / y_span               # legend height as axes fraction (0.5 data units)
        # Centre the legend at data y=-1 (middle of the 2-row leg_rows area).
        # axes_frac(y) = (y + leg_rows) / y_span
        ly = (leg_rows - 1.25) / y_span

        for cx_data, pal, low_lbl, high_lbl in [
            (acc_cx,                   acc_pal,  'Low/Poor', 'High/Good'),
            ((stab_cx + scale_cx) / 2, stab_pal, 'Low/Poor', 'High/Good'),
        ]:
            cx_frac = cx_data / total_width
            legend_ax = ax.inset_axes([cx_frac - lw / 2, ly, lw, lh])
            legend_ax.imshow(
                np.arange(len(pal)).reshape(1, len(pal)),
                cmap=mpl.colors.ListedColormap(list(pal)),
                interpolation='nearest', aspect='auto',
            )
            legend_ax.yaxis.set_ticks_position('none')
            legend_ax.xaxis.set_ticks_position('none')
            legend_ax.set_yticklabels([])
            legend_ax.set_xticks([0.5, len(pal) - 2])
            legend_ax.set_xticklabels([low_lbl, high_lbl], fontsize=9)
            for spine in legend_ax.spines.values():
                spine.set_visible(False)

        # ── Save ──────────────────────────────────────────────────────────────
        output_dir.mkdir(parents=True, exist_ok=True)
        stem = output_dir / 'Figure6'
        plt.savefig(stem.with_suffix('.pdf'), bbox_inches='tight')
        plt.savefig(stem.with_suffix('.png'), bbox_inches='tight')
        plt.close(fig)
        print(f"Saved Figure 6 to {stem}.pdf and .png")
