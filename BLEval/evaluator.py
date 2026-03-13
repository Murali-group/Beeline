from abc import ABC, abstractmethod
from pathlib import Path

import pandas as pd

from BLEval.data import EvaluationData


class Evaluator(ABC):
    """
    Abstract base class for BEELINE evaluation methods.

    Each subclass implements __call__ to compute a specific evaluation metric
    over an EvaluationData object and write results to disk. Output is written
    to each DatasetGroup's dataset_path so results are co-located with the
    predicted networks they describe.
    """

    @abstractmethod
    def __call__(self, evaluation_data: EvaluationData) -> None:
        """
        Compute the evaluation metric and write results to disk.

        Parameters
        ----------
        evaluation_data : EvaluationData
            Loaded predicted networks organised by dataset and run.

        Returns
        -------
        None
        """
        ...

    def _load_ground_truth(self, gt_path: Path) -> pd.DataFrame:
        """
        Load the ground truth network and return it as a DataFrame.

        Parameters
        ----------
        gt_path : Path
            Path to a CSV file with at least columns Gene1, Gene2.

        Returns
        -------
        pd.DataFrame
            Raw ground truth DataFrame with header row preserved.
        """
        if not isinstance(gt_path, Path):
            raise TypeError(f"gt_path must be Path, got {type(gt_path)}")

        return pd.read_csv(gt_path, header=0)

