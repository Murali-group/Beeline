from pathlib import Path
from typing import Dict, Iterator, List

import pandas as pd


class RunResult:
    """
    Predicted networks and ground truth reference for a single run.

    ranked_edges maps algorithm name to its ranked edge list DataFrame
    (columns: Gene1, Gene2, EdgeWeight). Algorithms whose rankedEdges.csv
    is missing are omitted with a warning rather than raising an error.
    """

    def __init__(
        self,
        run_id: str,
        ranked_edges: Dict[str, pd.DataFrame],
        ground_truth_path: Path,
        run_path: Path,
    ) -> None:
        # run_id : str  — matches run_id from the config 'runs' list
        # ranked_edges : dict[str, pd.DataFrame]  — keyed by algorithm_id
        # ground_truth_path : Path  — path to the ground truth edge list file
        # run_path : Path  — output directory for this run (output_dir/dataset_dir/dataset_id/run_id)
        if not isinstance(run_id, str):
            raise TypeError(f"run_id must be str, got {type(run_id)}")
        if not isinstance(ranked_edges, dict):
            raise TypeError(f"ranked_edges must be dict, got {type(ranked_edges)}")
        if not isinstance(ground_truth_path, Path):
            raise TypeError(f"ground_truth_path must be Path, got {type(ground_truth_path)}")
        if not isinstance(run_path, Path):
            raise TypeError(f"run_path must be Path, got {type(run_path)}")

        self.run_id: str = run_id
        self.ranked_edges: Dict[str, pd.DataFrame] = ranked_edges
        self.ground_truth_path: Path = ground_truth_path
        self.run_path: Path = run_path


class DatasetGroup:
    """
    All runs belonging to one dataset entry in the config.

    Runs in a DatasetGroup share the same ground truth network and represent
    multiple perturbations or noise realisations of the same biological system.
    Evaluation metrics (e.g. Jaccard, Spearman) are aggregated across runs
    within a group.
    """

    def __init__(self, dataset_id: str, runs: List[RunResult], dataset_path: Path) -> None:
        # dataset_id : str    — matches dataset_id from the config 'datasets' list
        # runs : list[RunResult]  — one entry per run in the config 'runs' list
        # dataset_path : Path — output directory for this dataset (output_dir/dataset_dir/dataset_id)
        if not isinstance(dataset_id, str):
            raise TypeError(f"dataset_id must be str, got {type(dataset_id)}")
        if not isinstance(runs, list) or not all(isinstance(r, RunResult) for r in runs):
            raise TypeError("runs must be a list of RunResult objects")
        if not isinstance(dataset_path, Path):
            raise TypeError(f"dataset_path must be Path, got {type(dataset_path)}")

        self.dataset_id: str = dataset_id
        self.runs: List[RunResult] = runs
        self.dataset_path: Path = dataset_path

    def __iter__(self) -> Iterator[RunResult]:
        return iter(self.runs)


class EvaluationData:
    """
    Loads and organises predicted networks from the output directory,
    mirroring the hierarchical structure of the 'datasets' section of config.yaml.

    Top level: datasets (DatasetGroup), each grouping multiple runs that share
    the same ground truth. Within each run, ranked edge lists are keyed by
    algorithm name. Algorithms with missing output files are skipped.
    """

    def __init__(self, config: dict, root: Path = None) -> None:
        # config : dict   — parsed YAML configuration
        # root   : Path   — working directory; defaults to cwd if None
        if not isinstance(config, dict):
            raise TypeError(f"config must be dict, got {type(config)}")

        if root is None:
            root = Path.cwd()
        if not isinstance(root, Path):
            raise TypeError(f"root must be Path, got {type(root)}")

        self.datasets: List[DatasetGroup] = self._load(config, root)

    def _load(self, config: dict, root: Path) -> List[DatasetGroup]:
        """
        Parse config and load rankedEdges DataFrames from disk.

        Parameters
        ----------
        config : dict
            Parsed YAML configuration dictionary.
        root : Path
            Root directory from which relative paths are resolved.

        Returns
        -------
        list[DatasetGroup]
            One DatasetGroup per enabled dataset entry in the config.
        """
        input_settings  = config['input_settings']
        output_settings = config['output_settings']

        input_dir   = root / input_settings['input_dir']
        output_dir  = root / output_settings['output_dir']
        dataset_dir = input_settings.get('dataset_dir', '')

        # Collect algorithm IDs that are enabled in the config.
        # should_run values may be bare booleans or single-element lists.
        algos: List[str] = [
            a['algorithm_id']
            for a in input_settings.get('algorithms', [])
            if (lambda v: v[0] if isinstance(v, list) else v)(
                a.get('should_run', [False])
            )
        ]

        groups: List[DatasetGroup] = []
        for ds in input_settings.get('datasets', []):
            should_run = ds.get('should_run', [True])
            if not (should_run[0] if isinstance(should_run, list) else should_run):
                continue

            dataset_id: str = ds['dataset_id']
            # outputs / dataset_dir / dataset_id
            dataset_path: Path = output_dir / dataset_dir / dataset_id
            runs: List[RunResult] = []

            for run in ds.get('runs', []):
                run_id: str = run['run_id']
                gt_filename: str = ds.get('groundTruthNetwork', 'GroundTruthNetwork.csv')

                # Ground truth: inputs / dataset_dir / dataset_id / filename
                gt_path: Path = input_dir / dataset_dir / dataset_id / gt_filename

                # outputs / dataset_dir / dataset_id / run_id
                run_path: Path = output_dir / dataset_dir / dataset_id / run_id

                # Ranked edges: run_path / algo / rankedEdges.csv
                ranked_edges: Dict[str, pd.DataFrame] = {}
                for algo in algos:
                    edges_path = run_path / algo / 'rankedEdges.csv'
                    if edges_path.exists():
                        ranked_edges[algo] = pd.read_csv(edges_path, sep='\t', header=0)
                    else:
                        print(f"Warning: {edges_path} not found, skipping.")

                runs.append(RunResult(run_id, ranked_edges, gt_path, run_path))

            groups.append(DatasetGroup(dataset_id, runs, dataset_path))

        return groups

    def __iter__(self) -> Iterator[DatasetGroup]:
        return iter(self.datasets)
