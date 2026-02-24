import argparse
import os
from pathlib import Path
import yaml
from tqdm import tqdm

from BLRun.genie3Runner import GENIE3Runner
from BLRun.grnboost2Runner import GRNBoost2Runner
from BLRun.grisliRunner import GRISLIRunner
from BLRun.grnvbemRunner import GRNVBEMRunner
from BLRun.jump3Runner import JUMP3Runner
from BLRun.leapRunner import LEAPRunner
from BLRun.pidcRunner import PIDCRunner
from BLRun.ppcorRunner import PPCORRunner
from BLRun.scodeRunner import SCODERunner
from BLRun.scribeRunner import SCRIBERunner
from BLRun.scsglRunner import SCSGLRunner
from BLRun.sinceritiesRunner import SINCERITIESRunner
from BLRun.singeRunner import SINGERunner

RUNNERS = {
    'GENIE3':       GENIE3Runner,
    'GRNBOOST2':    GRNBoost2Runner,
    'GRISLI':       GRISLIRunner,
    'GRNVBEM':      GRNVBEMRunner,
    'JUMP3':        JUMP3Runner,
    'LEAP':         LEAPRunner,
    'PIDC':         PIDCRunner,
    'PPCOR':        PPCORRunner,
    'SCODE':        SCODERunner,
    'SCRIBE':       SCRIBERunner,
    'SCSGL':        SCSGLRunner,
    'SINCERITIES':  SINCERITIESRunner,
    'SINGE':        SINGERunner,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description='BLRunner: Run GRN inference algorithms using BEELINE.'
    )
    parser.add_argument(
        '-c', '--config',
        type=str,
        required=True,
        help='Path to the configuration file used to run the inference algorithms.'
    )
    return parser.parse_args()


def load_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def get_datasets(input_settings):
    """
    Return a flat list of dataset dicts from input_settings.

    Each dataset entry in the config may have a 'variants' field (dict or list
    of dicts). Each variant becomes a separate runnable dataset. The dataset
    'dataset_id' is appended to dataset_dir so the runner resolves the correct path.

    If 'datasets' is absent, returns an empty list (caller handles auto-discovery).
    """
    if 'datasets' not in input_settings:
        return []

    datasets = []
    dataset_dir = input_settings.get('dataset_dir', '')

    for ds in input_settings['datasets']:
        if not ds.get('should_run', [True])[0]:
            continue

        if 'runs' not in ds:
            raise KeyError(f"Dataset '{ds['dataset_id']}' is missing required 'runs' field.")
        runs = ds['runs']
        # runs may be a single dict or a list of dicts
        if isinstance(runs, dict):
            runs = [runs]

        # Append dataset_id to dataset_dir so runner paths resolve correctly
        variant_dataset_dir = os.path.join(dataset_dir, ds['dataset_id']) if dataset_dir else ds['dataset_id']

        for run in runs:
            datasets.append({
                'dataset_dir': variant_dataset_dir,
                'dataset_id': run['run_id'],
                'exprData': run.get('exprData', 'ExpressionData.csv'),
                'pseudoTimeData': run.get('pseudoTimeData', 'PseudoTime.csv'),
                'groundTruthNetwork': ds.get('groundTruthNetwork', 'GroundTruthNetwork.csv'),
            })

    return datasets


def build_runner(algo_name, params, dataset, input_settings, output_settings):
    if algo_name not in RUNNERS:
        raise ValueError(f"Unknown algorithm '{algo_name}'. Available: {list(RUNNERS)}")

    runner_config = {
        'input': {
            'input_dir': input_settings['input_dir'],
            'dataset_dir': dataset['dataset_dir'],
        },
        'dataset': {
            'dataset_id':          dataset['dataset_id'],
            'exprData':            dataset['exprData'],
            'pseudoTimeData':      dataset['pseudoTimeData'],
            'groundTruthNetwork':  dataset['groundTruthNetwork'],
        },
        'output_settings': {
            'output_dir': output_settings['output_dir'],
            'run_id':     output_settings.get('run_id', ''),
        },
        'algo_name': algo_name,
        'params': params,
    }

    return RUNNERS[algo_name](Path.cwd(), runner_config)


def build_runners(config):
    input_settings  = config['input_settings']
    output_settings = config['output_settings']
    datasets        = get_datasets(input_settings)
    algorithms      = input_settings.get('algorithms', [])

    runners = []
    for dataset in datasets:
        for algo in algorithms:
            if not algo.get('should_run', [False])[0]:
                continue
            params = algo.get('params', {})
            runners.append(build_runner(algo['algorithm_id'], params, dataset, input_settings, output_settings))
    return runners


def get_working_dirs(config):
    """
    Compute expected working_dir paths from config without constructing Runner objects.

    Mirrors the path logic in Runner.__init__ so the overwrite check can run
    before runners are built (and their __init__ erases the directories).

    Parameters
    ----------
    config : dict
        Parsed YAML configuration dictionary.

    Returns
    -------
    list of Path
        One working_dir path per enabled dataset/algorithm combination.
    """
    root            = Path.cwd()
    input_settings  = config['input_settings']
    output_settings = config['output_settings']
    run_id          = output_settings.get('run_id', '')
    output_dir      = Path(output_settings['output_dir'])
    datasets        = get_datasets(input_settings)
    algorithms      = input_settings.get('algorithms', [])

    paths = []
    for dataset in datasets:
        for algo in algorithms:
            if not algo.get('should_run', [False])[0]:
                continue
            base_output = output_dir if output_dir.is_absolute() else root / output_dir
            if run_id:
                base_output = base_output / run_id
            if dataset['dataset_dir']:
                base_output = base_output / dataset['dataset_dir']
            base_output = base_output / dataset['dataset_id'] / algo['algorithm_id']
            paths.append(base_output / 'working_dir')

    return paths


def warn_if_populated(working_dirs):
    """
    Warn and prompt the user if any working directory already contains files.

    Parameters
    ----------
    working_dirs : list of Path
        Working directory paths to check for existing content.

    Returns
    -------
    bool
        True if the user confirms they want to proceed, False otherwise.
    """
    n_populated = sum(
        1 for p in working_dirs
        if p.exists() and any(p.iterdir())
    )
    if n_populated == 0:
        return True

    answer = input(
        f"Warning: {n_populated} working director{'y' if n_populated == 1 else 'ies'} "
        f"already exist and will be overwritten. Proceed? [y/n]: "
    ).strip().lower()
    return answer == 'y'


def main():
    args = parse_args()
    config = load_config(args.config)

    if not warn_if_populated(get_working_dirs(config)):
        print("Aborted.")
        return

    runners = build_runners(config)

    for runner in tqdm(runners):
        tqdm.write(runner.running_message)
        runner.generateInputs()
        runner.run()
        runner.parseOutput()


if __name__ == '__main__':
    main()
