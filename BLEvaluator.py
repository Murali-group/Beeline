import argparse

import yaml

from BLEval.AUPRC import AUPRC
from BLEval.AUROC import AUROC
from BLEval.data import EvaluationData
from BLEval.EarlyPrecision import EarlyPrecision
from BLEval.SignedEarlyPrecision import SignedEarlyPrecision


def parse_args():
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('-c', '--config', default='config.yaml',
        help="Configuration file containing list of datasets "
              "algorithms and output specifications.\n")

    parser.add_argument('-a', '--auc', action="store_true", default=False,
        help="Compute median of areas under Precision-Recall and ROC curves.\n")

    parser.add_argument('-j', '--jaccard', action="store_true", default=False,
        help="Compute median Jaccard index of predicted top-k networks "
             "for each algorithm for a given set of datasets generated "
             "from the same ground truth network.\n")

    parser.add_argument('-r', '--spearman', action="store_true", default=False,
        help="Compute median Spearman Corr. of predicted edges "
             "for each algorithm  for a given set of datasets generated "
             " from the same ground truth network.\n")

    parser.add_argument('-t', '--time', action="store_true", default=False,
        help="Analyze time taken by each algorithm for a.\n")

    parser.add_argument('-e', '--epr', action="store_true", default=False,
        help="Compute median early precision.")

    parser.add_argument('-s', '--sepr', action="store_true", default=False,
        help="Analyze median (signed) early precision for activation and inhibitory edges.")

    parser.add_argument('-m', '--motifs', action="store_true", default=False,
        help="Compute network motifs in the predicted top-k networks.")

    parser.add_argument('-p', '--paths', action="store_true", default=False,
        help="Compute path length statistics on the predicted top-k networks.")

    parser.add_argument('-b', '--borda', action="store_true", default=False,
        help="Compute edge ranked list using the various Borda aggregatio methods.")

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
    config = load_config(args.config)
    evaluation_data = EvaluationData(config)

    if args.auc:
        AUPRC()(evaluation_data)
        AUROC()(evaluation_data)

    if args.epr:
        EarlyPrecision()(evaluation_data)

    if args.sepr:
        SignedEarlyPrecision()(evaluation_data)


if __name__ == '__main__':
    main()
