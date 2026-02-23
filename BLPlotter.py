import argparse
import math
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import pandas as pd
import yaml

from BLPlot.plotter import box_plot, iter_datasets, load_metric
from BLPlot.PlotAUPRC import PlotAUPRC
from BLPlot.PlotAUROC import PlotAUROC
from BLPlot.PlotSummaryHeatmap import PlotSummaryHeatmap
from BLPlot.PlotEPRHeatmap import PlotEPRHeatmap


def parse_args():
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Plot BEELINE evaluation results.')

    parser.add_argument('-c', '--config', required=True,
        help="Configuration file containing list of datasets and algorithms. "
             "The same file used with BLEvaluator.py may be used here.\n")

    parser.add_argument('-o', '--output', default='./',
        help="Output directory for generated plots.\n")

    parser.add_argument('-a', '--auprc', action='store_true', default=False,
        help="Produce a box plot of AUPRC values for the evaluated algorithms.\n")

    parser.add_argument('-r', '--auroc', action='store_true', default=False,
        help="Produce a box plot of AUROC values for the evaluated algorithms.\n")

    parser.add_argument('-e', '--epr', action='store_true', default=False,
        help="Produce a box plot of early precision values for the evaluated algorithms.\n")

    parser.add_argument('-v', '--overview', action='store_true', default=False,
        help="Produce a combined plot of AUPRC and early precision ratios.\n")

    parser.add_argument('--summary', action='store_true', default=False,
        help="Produce a Figure-2-style heatmap of median AUPRC ratios and "
             "median Spearman stability per algorithm and dataset.\n")

    parser.add_argument('--epr-summary', action='store_true', default=False,
        help="Produce a Figure-4-style heatmap of median AUPRC ratio, EPR "
             "ratio, and signed EPR ratios per algorithm and dataset.\n")

    return parser.parse_args()


def load_config(config_path: str) -> dict:
    """
    Load and return the YAML configuration file.

    Parameters
    ----------
    config_path : str
        Path to the YAML configuration file.

    Returns
    -------
    dict
        Parsed configuration dictionary.
    """
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def _random_epr(gt_path: Path) -> float:
    """
    Compute the expected early precision of a random predictor for a dataset.

    The random baseline is k / (n*(n-1)) where k is the number of non-self-loop
    ground truth edges and n is the number of unique genes.

    Parameters
    ----------
    gt_path : Path
        Path to the ground truth edge list CSV (columns Gene1, Gene2).

    Returns
    -------
    float
        Expected early precision of a random predictor, or nan if undefined.
    """
    if not isinstance(gt_path, Path):
        raise TypeError(f"gt_path must be Path, got {type(gt_path)}")

    gt = pd.read_csv(gt_path, header=0)
    # Exclude self-loops to match the EPR evaluator's definition of k
    gt = gt[gt['Gene1'] != gt['Gene2']]
    k = len(gt)

    genes = set(gt['Gene1']).union(set(gt['Gene2']))
    n = len(genes)
    n_possible = n * (n - 1)

    if n_possible == 0 or k == 0:
        return float('nan')

    return k / n_possible


def _load_epr_ratios(
    config: dict,
    root: Path,
) -> Dict[str, List[float]]:
    """
    Load early precision values and normalize each by the random predictor baseline.

    EPR ratio = early_precision / (k / (n*(n-1))), where k is the number of
    non-self-loop ground truth edges and n is the number of unique genes in the
    dataset. Datasets with missing EarlyPrecision.csv or missing ground truth are
    skipped with a warning.

    Parameters
    ----------
    config : dict
        Parsed YAML configuration.
    root : Path
        Working directory from which config paths are resolved.

    Returns
    -------
    dict[str, list[float]]
        Algorithm name -> list of EPR ratio values across all datasets/runs.
    """
    if not isinstance(config, dict):
        raise TypeError(f"config must be dict, got {type(config)}")
    if not isinstance(root, Path):
        raise TypeError(f"root must be Path, got {type(root)}")

    all_values: Dict[str, List[float]] = {}

    for _, dataset_path, gt_path in iter_datasets(config, root):
        csv_path = dataset_path / 'EarlyPrecision.csv'
        if not csv_path.exists():
            print(f"Warning: {csv_path} not found, skipping.")
            continue
        if not gt_path.exists():
            print(f"Warning: ground truth {gt_path} not found, skipping.")
            continue

        baseline = _random_epr(gt_path)
        if math.isnan(baseline) or baseline == 0.0:
            print(f"Warning: random EPR is undefined for {gt_path}, skipping.")
            continue

        df = pd.read_csv(csv_path, index_col=0)
        for algo in df.index:
            ratios = [
                v / baseline
                for v in df.loc[algo].tolist()
                if not math.isnan(v)
            ]
            all_values.setdefault(str(algo), []).extend(ratios)

    return all_values


def _overview_plot(
    auprc_values: Dict[str, List[float]],
    epr_ratio_values: Dict[str, List[float]],
    output_path: Path,
) -> None:
    """
    Save a two-panel overview plot of AUPRC and EPR ratios side by side.

    Both panels share the same set of algorithms on the x-axis. Algorithms
    present in only one panel appear as empty boxes in the other.

    Parameters
    ----------
    auprc_values : dict[str, list[float]]
        Algorithm name -> list of AUPRC values.
    epr_ratio_values : dict[str, list[float]]
        Algorithm name -> list of EPR ratio values.
    output_path : Path
        Destination file path (PDF).
    """
    if not isinstance(auprc_values, dict):
        raise TypeError(f"auprc_values must be dict, got {type(auprc_values)}")
    if not isinstance(epr_ratio_values, dict):
        raise TypeError(f"epr_ratio_values must be dict, got {type(epr_ratio_values)}")
    if not isinstance(output_path, Path):
        raise TypeError(f"output_path must be Path, got {type(output_path)}")

    algos = sorted(set(auprc_values) | set(epr_ratio_values))
    if not algos:
        print("No data to plot for overview.")
        return

    fig, (ax_auprc, ax_epr) = plt.subplots(1, 2, figsize=(max(10, len(algos) * 1.6), 5))

    ax_auprc.boxplot([auprc_values.get(a, []) for a in algos], labels=algos)
    ax_auprc.set_title('AUPRC')
    ax_auprc.set_ylabel('AUPRC')
    ax_auprc.set_xlabel('Algorithm')
    plt.setp(ax_auprc.get_xticklabels(), rotation=45, ha='right')

    ax_epr.boxplot([epr_ratio_values.get(a, []) for a in algos], labels=algos)
    ax_epr.set_title('Early Precision Ratio')
    ax_epr.set_ylabel('EPR / Random EPR')
    ax_epr.set_xlabel('Algorithm')
    plt.setp(ax_epr.get_xticklabels(), rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close(fig)
    print(f"Saved overview plot to {output_path}")


def main():
    args   = parse_args()
    config = load_config(args.config)
    root   = Path.cwd()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.auprc:
        PlotAUPRC()(config, output_dir, root)

    if args.auroc:
        PlotAUROC()(config, output_dir, root)

    if args.epr:
        values = load_metric(config, 'EarlyPrecision.csv', root)
        box_plot(values, 'Early Precision', 'Early Precision', output_dir / 'EarlyPrecision.pdf')

    if args.summary:
        PlotSummaryHeatmap()(config, output_dir, root)

    if args.epr_summary:
        PlotEPRHeatmap()(config, output_dir, root)

    if args.overview:
        auprc_values     = load_metric(config, 'AUPRC.csv', root)
        epr_ratio_values = _load_epr_ratios(config, root)
        _overview_plot(auprc_values, epr_ratio_values, output_dir / 'Overview.pdf')


if __name__ == '__main__':
    main()
