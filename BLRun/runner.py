from abc import ABC, abstractmethod
from pathlib import Path

class Runner(ABC):
    """
    Abstract base_input class for BEELINE GRN inference algorithm runners.

    Subclasses must implement generateInputs, run, and parseOutput.
    Attributes set here reflect the fields accessed by runner implementations.
    """

    def __init__(self, root: Path, config: dict):
        """
        Parameters
        ----------
        root : Path
            Root path from which all subpaths are resolved.
        config : dict
            Merged configuration for a single dataset + algorithm run.
            Expected structure:
              input:
                input_dir:   <str>  input directory (absolute, or relative to root)
                dataset_dir: <str>  optional subdirectory (may be empty)
              dataset:
                name:                <str>  subdirectory name / dataset label
                exprData:            <str>  expression data filename
                pseudoTimeData:      <str>  pseudotime data filename
                groundTruthNetwork:  <str>  ground truth network filename
              output_settings:
                output_dir: <str>  output directory (absolute, or relative to root)
              algo_name: <str>  name of the algorithm (appended to output_dir)
              params: <dict>  algorithm-specific parameters
        """
        inp = config['input']
        ds  = config['dataset']

        input_dir_path  = Path(inp['input_dir'])
        output_dir_path = Path(config['output_settings']['output_dir'])

        base_input = input_dir_path if input_dir_path.is_absolute() else root / input_dir_path
        if inp.get('dataset_dir'):
            base_input = base_input / inp['dataset_dir']

        base_output = output_dir_path if output_dir_path.is_absolute() else root / output_dir_path
        if inp.get('dataset_dir'):
            base_output = base_output / inp['dataset_dir']

        self.input_dir   = base_input / ds['name']
        self.output_dir  = base_output / ds['name']
        self.working_dir = self.output_dir / config['algo_name']
        self.exprData           = ds.get('exprData',           'ExpressionData.csv')
        self.pseudoTimeData     = ds.get('pseudoTimeData',     'PseudoTime.csv')
        self.groundTruthNetwork = ds.get('groundTruthNetwork', 'GroundTruthNetwork.csv')
        self.params             = config.get('params', {})

    @abstractmethod
    def generateInputs(self):
        """Prepare algorithm-specific input files from the dataset."""

    @abstractmethod
    def run(self):
        """Execute the inference algorithm."""

    @abstractmethod
    def parseOutput(self):
        """Parse algorithm output into a ranked edge list (rankedEdges.csv)."""
