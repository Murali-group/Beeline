import argparse
from pathlib import Path

import yaml

from BLPlot.PlotAUPRC import PlotAUPRC
from BLPlot.PlotAUROC import PlotAUROC
from BLPlot.PlotSummaryHeatmap import PlotSummaryHeatmap
from BLPlot.PlotEPRHeatmap import PlotEPRHeatmap
from BLPlot.PlotEPR import PlotEPR
from BLPlot.PlotFigure5 import PlotFigure5
from BLPlot.PlotFigure6 import PlotFigure6


def parse_args():
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Plot BEELINE evaluation results.')

    parser.add_argument('-c', '--config', default=None,
        help="Configuration file containing list of datasets and algorithms. "
             "The same file used with BLEvaluator.py may be used here. "
             "Required for all plots except --figure-6.\n")

    parser.add_argument('-o', '--output', default='./',
        help="Output directory for generated plots.\n")

    parser.add_argument('-a', '--auprc', action='store_true', default=False,
        help="Produce per-dataset AUPRC plots (AUPRC.pdf). Datasets with a "
             "single run output a precision-recall curve; datasets with "
             "multiple runs output a box plot.\n")

    parser.add_argument('-r', '--auroc', action='store_true', default=False,
        help="Produce per-dataset AUROC plots (AUROC.pdf). Datasets with a "
             "single run output a ROC curve; datasets with multiple runs "
             "output a box plot.\n")

    parser.add_argument('-e', '--epr', action='store_true', default=False,
        help="Produce a box plot of early precision values for the evaluated "
             "algorithms (EPR.pdf).\n")

    parser.add_argument('--summary', action='store_true', default=False,
        help="Produce a Figure-2-style heatmap of median AUPRC ratios and "
             "median Spearman stability per algorithm and dataset (Summary.pdf).\n")

    parser.add_argument('--epr-summary', action='store_true', default=False,
        help="Produce a Figure-4-style heatmap of median AUPRC ratio, EPR "
             "ratio, and signed EPR ratios per algorithm and dataset "
             "(EPRSummary.pdf).\n")

    parser.add_argument('--figure-5', action='store_true', default=False,
        help="Produce a Figure-5-style table of network statistics and EPR "
             "values, pairing TFs + 500 genes (left half) and TFs + 1000 "
             "genes (right half) for each dataset (Figure5.pdf).\n")

    parser.add_argument('--figure-6', metavar='YAML', default=None,
        help="Produce a Figure-6-style algorithm comparison table from a "
             "standalone YAML config file (Figure6.pdf). Mutually exclusive "
             "with -c/--config.\n")

    parser.add_argument('--all', action='store_true', default=False,
        help="Run all plots.\n")

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


def main():
    args = parse_args()

    if args.figure_6 and args.config:
        raise SystemExit("error: --figure-6 and -c/--config are mutually exclusive.")

    if args.figure_6:
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        fig6_config = load_config(args.figure_6)
        PlotFigure6()(fig6_config, output_dir)
        return

    if not args.config:
        raise SystemExit("error: -c/--config is required unless --figure-6 is used.")

    config = load_config(args.config)
    root   = Path.cwd()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_all = args.all

    if args.auprc or run_all:
        PlotAUPRC()(config, output_dir, root)

    if args.auroc or run_all:
        PlotAUROC()(config, output_dir, root)

    if args.epr or run_all:
        PlotEPR()(config, output_dir, root)

    if args.summary or run_all:
        PlotSummaryHeatmap()(config, output_dir, root)

    if args.epr_summary or run_all:
        PlotEPRHeatmap()(config, output_dir, root)

    if args.figure_5 or run_all:
        PlotFigure5()(config, output_dir, root)


if __name__ == '__main__':
    main()
