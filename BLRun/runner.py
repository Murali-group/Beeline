from abc import ABC, abstractmethod
from pathlib import Path
import subprocess
import pandas as pd

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
              dataset:
                dataset_id:          <str>  dataset group label (path segment under input_dir)
                run_id:              <str>  run label (path segment under dataset_id)
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
        # experiment_id : str — optional; when set, an experiment_id segment is
        # inserted between output_dir and the dataset path so multiple experiment
        # runs can coexist under the same base output directory.
        experiment_id_prefix = config['output_settings'].get('experiment_id', '')

        base_input = input_dir_path if input_dir_path.is_absolute() else root / input_dir_path

        base_output = output_dir_path if output_dir_path.is_absolute() else root / output_dir_path
        if experiment_id_prefix:
            base_output = base_output / experiment_id_prefix
        base_output = base_output / ds['dataset_id'] / ds['run_id'] / config['algo_name']

        base_input = base_input.resolve()
        base_output = base_output.resolve()

        # input_dir: run-level input directory (expression data, pseudo-time).
        self.input_dir  = base_input / ds['dataset_id'] / ds['run_id']
        self.output_dir = base_output
        self.working_dir = base_output / "working_dir"

        # Erase working directory so stale inputs from prior runs are not reused.
        if self.working_dir.exists():
            for item in sorted(self.working_dir.rglob('*'), reverse=True):
                item.unlink() if (item.is_file() or item.is_symlink()) else item.rmdir()
            self.working_dir.rmdir()

        # Pre-create output_dir and working_dir so docker cannot claim them as root.
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.working_dir.mkdir(parents=True, exist_ok=True)
        
        # Precompute progress message for CLI output.
        self.running_message = (
            f"Running {config['algo_name']} | dataset: {ds['dataset_id']} | run: {ds['run_id']}"
        )

        self.exprData           = ds.get('exprData',           'ExpressionData.csv')
        self.pseudoTimeData     = ds.get('pseudoTimeData',     'PseudoTime.csv')
        self.groundTruthNetwork = ds.get('groundTruthNetwork', 'GroundTruthNetwork.csv')
        # ground_truth_file: full path to the dataset-level ground truth CSV.
        self.ground_truth_file  = base_input / ds['dataset_id'] / self.groundTruthNetwork
        
        # image: Docker image name used to run this algorithm (e.g. "grnbeeline/genie3:base").
        # Mandatory — every algorithm entry in the config must supply this field.
        if 'image' not in config or not config['image']:
            raise ValueError("Algorithm config must include a non-empty 'image' field.")
        if not isinstance(config['image'], str):
            raise TypeError(f"'image' must be a str, got {type(config['image'])}")
        self.image = config['image']

        # Unwrap single-element lists so runners receive scalar values.
        # YAML config files commonly wrap param values in brackets
        # (e.g. `pVal: [0.01]`), which YAML parses as a list.
        raw_params = config.get('params', {})
        self.params = {
            k: (v[0] if isinstance(v, list) and len(v) == 1 else v)
            for k, v in raw_params.items()
        }

    @abstractmethod
    def generateInputs(self):
        """Prepare algorithm-specific input files from the dataset."""

    @abstractmethod
    def run(self):
        """Execute the inference algorithm."""

    @abstractmethod
    def parseOutput(self) -> None:
        """
        Parse raw algorithm output and write a ranked edge list to disk.

        Implementations should build a DataFrame with columns Gene1, Gene2,
        EdgeWeight and pass it to self._write_ranked_edges(). Returns early
        without writing if the expected output file is missing.
        """

    def _run_docker(self, cmd: str, append: bool = False) -> None:
        """
        Execute a shell command and write combined stdout/stderr to output.txt.

        Parameters
        ----------
        cmd : str
            Shell command to execute (passed to the shell verbatim).
        append : bool
            If True, append to an existing output.txt. Use for runners that
            invoke docker in a loop so all container output is collected in
            one file. Defaults to False (overwrite).
        """
        if not isinstance(cmd, str):
            raise TypeError(f"cmd must be str, got {type(cmd)}")
        if not isinstance(append, bool):
            raise TypeError(f"append must be bool, got {type(append)}")

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        mode = 'a' if append else 'w'
        with open(self.output_dir / 'output.txt', mode) as f:
            f.write(result.stdout)
            f.write(result.stderr)

        if result.returncode != 0:
            raise RuntimeError(
                f"Docker command failed (exit {result.returncode}). "
                f"See {self.output_dir / 'output.txt'} for details."
            )

    def _write_ranked_edges(self, df: pd.DataFrame) -> None:
        """
        Write a ranked edge list to self.output_dir/rankedEdges.csv.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with columns Gene1, Gene2, EdgeWeight.

        Returns
        -------
        None

        Raises
        ------
        FileNotFoundError
            If self.output_dir does not exist at the time of writing.
        TypeError
            If df is not a pd.DataFrame.
        """
        if not isinstance(df, pd.DataFrame):
            raise TypeError(f"df must be pd.DataFrame, got {type(df)}")
        if not self.output_dir.is_dir():
            raise FileNotFoundError(
                f"Output directory does not exist: {self.output_dir}")
        df.to_csv(self.output_dir / 'rankedEdges.csv', sep='\t', index=False)
