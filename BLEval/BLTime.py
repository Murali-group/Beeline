import re
from pathlib import Path
from typing import Dict, List

import pandas as pd

from BLEval.evaluator import Evaluator
from BLEval.data import EvaluationData


def _parse_cpu_time(time_file: Path) -> float:
    """
    Parse CPU time from a single GNU ``time -v`` output file.

    CPU time is defined as the sum of user time and system time reported by
    ``time -v``. Both fields are required; if either is missing the file is
    considered malformed and 0.0 is returned with a warning.

    Parameters
    ----------
    time_file : Path
        Path to a file produced by ``time -v -o <file> <command>``.

    Returns
    -------
    float
        CPU time in seconds (user + system).
    """
    if not isinstance(time_file, Path):
        raise TypeError(f"time_file must be Path, got {type(time_file)}")

    user_time: float = 0.0
    sys_time: float = 0.0
    found_user = False
    found_sys = False

    with open(time_file) as fh:
        for line in fh:
            m = re.search(r'User time \(seconds\):\s+([\d.]+)', line)
            if m:
                user_time = float(m.group(1))
                found_user = True
                continue
            m = re.search(r'System time \(seconds\):\s+([\d.]+)', line)
            if m:
                sys_time = float(m.group(1))
                found_sys = True

    if not (found_user and found_sys):
        print(f"Warning: could not parse CPU time from {time_file}, treating as 0.0.")
        return 0.0

    return user_time + sys_time


def _total_cpu_time(working_dir: Path) -> float:
    """
    Sum CPU times across all ``time*.txt`` files in working_dir.

    Multi-trajectory algorithms produce one timing file per trajectory
    (``time0.txt``, ``time1.txt``, …); this function accumulates them so the
    return value represents the total CPU time for the full algorithm run.
    Returns float('nan') when no timing files are present.

    Parameters
    ----------
    working_dir : Path
        Algorithm working directory (``run_path / algo / working_dir``).

    Returns
    -------
    float
        Total CPU time in seconds, or nan if no timing files exist.
    """
    if not isinstance(working_dir, Path):
        raise TypeError(f"working_dir must be Path, got {type(working_dir)}")

    time_files: List[Path] = sorted(working_dir.glob('time*.txt'))

    if not time_files:
        return float('nan')

    return sum(_parse_cpu_time(f) for f in time_files)


class BLTime(Evaluator):
    """
    Evaluator that reports the CPU time consumed by each algorithm.

    Timing files (``time*.txt``) are produced by the ``time -v`` shell
    utility and written to each algorithm's ``working_dir`` during the run
    phase. CPU time is defined as user time + system time; multi-trajectory
    algorithms may produce multiple files whose values are summed.

    For each DatasetGroup, writes ``time.csv`` to ``dataset_path``. Rows are
    algorithms and columns are run_ids. Missing timing files produce nan
    entries. A dictionary mapping algorithm name to CPU time (in seconds) is
    returned for each dataset group.
    """

    def __call__(self, evaluation_data: EvaluationData) -> None:
        """
        Compute CPU time per algorithm per run and write results to
        ``dataset_path/time.csv``.

        For each run, every algorithm's ``working_dir`` is searched for
        ``time*.txt`` files. Their CPU times are summed to give one value per
        algorithm per run. Results are written as a CSV with algorithms as
        rows and run_ids as columns.

        Parameters
        ----------
        evaluation_data : EvaluationData
            Loaded predicted networks organised by dataset and run.

        Returns
        -------
        None
        """
        if not isinstance(evaluation_data, EvaluationData):
            raise TypeError(
                f"evaluation_data must be EvaluationData, got {type(evaluation_data)}"
            )

        for dataset_group in evaluation_data:
            # results[algo][run_id] = total CPU time in seconds
            results: Dict[str, Dict[str, float]] = {}

            for run in dataset_group:
                # Collect all algo names visible in this run's ranked_edges,
                # plus any algos that only have timing dirs but no ranked output.
                algo_dirs = {
                    d.name
                    for d in run.run_path.iterdir()
                    if d.is_dir()
                } if run.run_path.exists() else set()

                algos = set(run.ranked_edges.keys()) | algo_dirs

                for algo in algos:
                    working_dir = run.run_path / algo / 'working_dir'
                    cpu_time = _total_cpu_time(working_dir)
                    results.setdefault(algo, {})[run.run_id] = cpu_time

            if not results:
                continue

            # Build output DataFrame: rows = algorithms, columns = run_ids
            out_df = pd.DataFrame(results).T
            out_df.index.name = 'Algorithm'

            dataset_group.dataset_path.mkdir(parents=True, exist_ok=True)
            out_path = dataset_group.dataset_path / 'time.csv'
            out_df.to_csv(out_path)
            print(f"Wrote timing results to {out_path}")
