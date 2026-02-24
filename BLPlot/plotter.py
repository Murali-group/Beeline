import math
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Iterator, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def get_algo_ids(config: dict) -> List[str]:
    """
    Return the list of enabled algorithm IDs from the config.

    Parameters
    ----------
    config : dict
        Parsed YAML configuration.

    Returns
    -------
    list[str]
        Algorithm IDs whose should_run flag is truthy.
    """
    if not isinstance(config, dict):
        raise TypeError(f"config must be dict, got {type(config)}")

    algos = []
    for a in config['input_settings'].get('algorithms', []):
        flag = a.get('should_run', [False])
        if flag[0] if isinstance(flag, list) else flag:
            algos.append(a['algorithm_id'])
    return algos


def iter_datasets_with_runs(
    config: dict,
    root: Path,
) -> Iterator[Tuple[str, Path, Path, list]]:
    """
    Yield (dataset_id, dataset_path, gt_path, runs) for each enabled dataset.

    Extends iter_datasets by also yielding the raw list of run dicts from the
    config, so callers can determine run count and resolve per-run paths via
    dataset_path / run['run_id'].

    Parameters
    ----------
    config : dict
        Parsed YAML configuration.
    root : Path
        Working directory from which config paths are resolved.

    Yields
    ------
    tuple of (str, Path, Path, list)
        dataset_id, output dataset directory, ground truth file path,
        list of run dicts from config (each has at least 'run_id').
    """
    if not isinstance(config, dict):
        raise TypeError(f"config must be dict, got {type(config)}")
    if not isinstance(root, Path):
        raise TypeError(f"root must be Path, got {type(root)}")

    input_settings  = config['input_settings']
    output_settings = config['output_settings']

    input_dir  = root / input_settings['input_dir']
    output_dir = root / output_settings['output_dir']
    # run_id : str — optional; when set, a run_id segment is inserted
    # between output_dir and the dataset path.
    run_id = output_settings.get('run_id', '')
    if run_id:
        output_dir = output_dir / run_id
    dataset_dir = input_settings.get('dataset_dir', '')

    for ds in input_settings.get('datasets', []):
        should_run = ds.get('should_run', [True])
        if not (should_run[0] if isinstance(should_run, list) else should_run):
            continue

        dataset_id  = ds['dataset_id']
        gt_filename = ds.get('groundTruthNetwork', 'GroundTruthNetwork.csv')

        dataset_path = output_dir / dataset_dir / dataset_id
        gt_path      = input_dir  / dataset_dir / dataset_id / gt_filename
        runs         = ds.get('runs', [])

        yield dataset_id, dataset_path, gt_path, runs


def iter_datasets(
    config: dict,
    root: Path,
) -> Iterator[Tuple[str, Path, Path]]:
    """
    Yield (dataset_id, dataset_path, gt_path) for each enabled dataset in config.

    Parameters
    ----------
    config : dict
        Parsed YAML configuration.
    root : Path
        Working directory from which config paths are resolved.

    Yields
    ------
    tuple of (str, Path, Path)
        dataset_id, output dataset directory, ground truth file path.
    """
    if not isinstance(config, dict):
        raise TypeError(f"config must be dict, got {type(config)}")
    if not isinstance(root, Path):
        raise TypeError(f"root must be Path, got {type(root)}")

    input_settings  = config['input_settings']
    output_settings = config['output_settings']

    input_dir  = root / input_settings['input_dir']
    output_dir = root / output_settings['output_dir']
    # run_id : str — optional; when set, a run_id segment is inserted
    # between output_dir and the dataset path.
    run_id = output_settings.get('run_id', '')
    if run_id:
        output_dir = output_dir / run_id
    dataset_dir = input_settings.get('dataset_dir', '')

    for ds in input_settings.get('datasets', []):
        should_run = ds.get('should_run', [True])
        if not (should_run[0] if isinstance(should_run, list) else should_run):
            continue

        dataset_id  = ds['dataset_id']
        gt_filename = ds.get('groundTruthNetwork', 'GroundTruthNetwork.csv')

        dataset_path = output_dir / dataset_dir / dataset_id
        gt_path      = input_dir  / dataset_dir / dataset_id / gt_filename

        yield dataset_id, dataset_path, gt_path


def load_metric(
    config: dict,
    metric_csv: str,
    root: Path,
) -> Dict[str, List[float]]:
    """
    Aggregate per-algorithm metric values from pre-computed CSV files.

    Reads {dataset_path}/{metric_csv} for each enabled dataset. The CSV is
    expected to have rows = algorithms (index column 'Algorithm') and columns
    = run_ids. Missing files are skipped with a warning.

    Parameters
    ----------
    config : dict
        Parsed YAML configuration.
    metric_csv : str
        Filename of the metric CSV (e.g. 'AUPRC.csv').
    root : Path
        Working directory from which config paths are resolved.

    Returns
    -------
    dict[str, list[float]]
        Algorithm name -> list of all non-NaN values across all datasets/runs.
    """
    if not isinstance(config, dict):
        raise TypeError(f"config must be dict, got {type(config)}")
    if not isinstance(metric_csv, str):
        raise TypeError(f"metric_csv must be str, got {type(metric_csv)}")
    if not isinstance(root, Path):
        raise TypeError(f"root must be Path, got {type(root)}")

    all_values: Dict[str, List[float]] = {}

    for _, dataset_path, _ in iter_datasets(config, root):
        csv_path = dataset_path / metric_csv
        if not csv_path.exists():
            print(f"Warning: {csv_path} not found, skipping.")
            continue

        df = pd.read_csv(csv_path, index_col=0)
        for algo in df.index:
            values = [v for v in df.loc[algo].tolist() if not math.isnan(v)]
            all_values.setdefault(str(algo), []).extend(values)

    return all_values


def random_classifier_baseline(gt_path: Path) -> float:
    """
    Compute the expected precision of a random predictor for a dataset.

    Equals k / (n*(n-1)) where k is the number of non-self-loop ground truth
    edges and n is the number of unique genes. This is the random baseline for
    both AUPRC and early precision (the PR curve of a random predictor is flat
    at height k / n_possible).

    Parameters
    ----------
    gt_path : Path
        Path to the ground truth edge list CSV (columns Gene1, Gene2).

    Returns
    -------
    float
        Random predictor baseline in (0, 1], or nan if undefined (empty network).
    """
    if not isinstance(gt_path, Path):
        raise TypeError(f"gt_path must be Path, got {type(gt_path)}")

    gt = pd.read_csv(gt_path, header=0)
    gt = gt[gt['Gene1'] != gt['Gene2']]
    k = len(gt)

    genes = set(gt['Gene1']).union(set(gt['Gene2']))
    n = len(genes)
    n_possible = n * (n - 1)

    if n_possible == 0 or k == 0:
        return float('nan')

    return k / n_possible


def load_dataset_metric(
    dataset_path: Path,
    metric_csv: str,
) -> Dict[str, List[float]]:
    """
    Load per-algorithm metric values from one dataset's pre-computed CSV.

    The CSV is expected to have rows = algorithms (index column) and columns
    = run_ids. Returns an empty dict (with a warning) if the file is missing.

    Parameters
    ----------
    dataset_path : Path
        Output directory for the dataset (contains the metric CSV).
    metric_csv : str
        Filename of the metric CSV (e.g. 'AUPRC.csv').

    Returns
    -------
    dict[str, list[float]]
        Algorithm name -> list of non-NaN run values for this dataset.
    """
    if not isinstance(dataset_path, Path):
        raise TypeError(f"dataset_path must be Path, got {type(dataset_path)}")
    if not isinstance(metric_csv, str):
        raise TypeError(f"metric_csv must be str, got {type(metric_csv)}")

    csv_path = dataset_path / metric_csv
    if not csv_path.exists():
        print(f"Warning: {csv_path} not found, skipping.")
        return {}

    df = pd.read_csv(csv_path, index_col=0)
    return {
        str(algo): [v for v in df.loc[algo].tolist() if not math.isnan(v)]
        for algo in df.index
    }


def make_box_figure(
    values: Dict[str, List[float]],
    title: str,
    ylabel: str,
    rand_value: float = None,
) -> 'plt.Figure | None':
    """
    Create and return a box plot figure without saving it.

    Renders a seaborn box plot with individual data points overlaid as a swarm
    plot. When rand_value is provided, a dashed grey reference line marks the
    expected performance of a random predictor. Returns None when values is
    empty so callers can skip adding an empty page to a multi-page PDF.

    Parameters
    ----------
    values : dict[str, list[float]]
        Algorithm name -> list of observed metric values.
    title : str
        Plot title.
    ylabel : str
        Y-axis label.
    rand_value : float or None
        If provided, a dashed grey horizontal line is drawn at this y value to
        indicate the random-predictor baseline.

    Returns
    -------
    plt.Figure or None
        The created figure, or None if values is empty.
    """
    if not isinstance(values, dict):
        raise TypeError(f"values must be dict, got {type(values)}")

    if not values:
        print(f"No data to plot for '{title}'.")
        return None

    plt.rcParams.update({'font.size': 14})

    algos = sorted(values.keys())

    # Build long-form DataFrame required by seaborn
    records = [(algo, v) for algo in algos for v in values[algo]]
    df = pd.DataFrame(records, columns=['Algorithm', 'Value'])

    fig, ax = plt.subplots(figsize=(max(6, len(algos) * 0.9), 5))

    # Dashed grey line marking the random-predictor baseline
    if rand_value is not None:
        ax.axhline(rand_value, color='gray', linestyle='--', linewidth=0.8)

    sns.boxplot(
        data=df, x='Algorithm', y='Value', order=algos,
        palette=sns.color_palette("Set1", n_colors=len(algos)),
        fliersize=0, ax=ax,
    )
    sns.swarmplot(
        data=df, x='Algorithm', y='Value', order=algos,
        alpha=0.5, color='k', ax=ax,
    )

    ax.set_ylim([0.0, 1.0])
    ax.set_title(title)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_xlabel('Algorithm', fontsize=18)
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    plt.tight_layout()
    return fig


def box_plot(
    values: Dict[str, List[float]],
    title: str,
    ylabel: str,
    output_path: Path,
) -> None:
    """
    Save a single-page box plot PDF.

    Convenience wrapper around make_box_figure for callers that produce one
    plot per file. The figure is closed after saving.

    Parameters
    ----------
    values : dict[str, list[float]]
        Algorithm name -> list of observed metric values.
    title : str
        Plot title.
    ylabel : str
        Y-axis label.
    output_path : Path
        Destination file path (PDF).
    """
    if not isinstance(output_path, Path):
        raise TypeError(f"output_path must be Path, got {type(output_path)}")

    fig = make_box_figure(values, title, ylabel)
    if fig is None:
        return
    fig.savefig(output_path)
    plt.close(fig)
    print(f"Saved plot to {output_path}")


class Plotter(ABC):
    """
    Abstract base class for BEELINE plot generators.

    Each subclass implements __call__ to read pre-computed evaluation CSVs and
    write one or more plot files to a caller-specified output directory.
    Shared loading and rendering helpers are provided as module-level functions
    in this module: iter_datasets, load_metric, box_plot.
    """

    @abstractmethod
    def __call__(self, config: dict, output_dir: Path, root: Path) -> None:
        """
        Generate plots from pre-computed evaluation CSVs.

        Parameters
        ----------
        config : dict
            Parsed YAML configuration.
        output_dir : Path
            Directory where plot files are written.
        root : Path
            Working directory from which config paths are resolved.

        Returns
        -------
        None
        """
        ...
