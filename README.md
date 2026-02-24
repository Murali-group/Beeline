# BEELINE: Benchmarking gEnE reguLatory network Inference from siNgle-cEll transcriptomic data

![Overview of BEELINE](docs/figs/overview-graphic.png)

BEELINE is a benchmarking framework for evaluating gene regulatory network (GRN) inference algorithms on single-cell RNA-seq data. It runs algorithms via Docker containers, evaluates their output against a ground truth network, and produces summary plots.

Full documentation: [https://murali-group.github.io/Beeline/](https://murali-group.github.io/Beeline/)

## Setup

**1. Install the conda environment**
```bash
bash setupAnacondaVENV.sh
```

**2. Pull algorithm Docker images**
```bash
bash initialize.sh
```
Pass the `--build` flag to build images locally instead of pulling from DockerHub.

**3. Activate the environment**
```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate BEELINE
```

---

## Usage

All three pipeline stages take a YAML configuration file via `-c/--config`.

### 1. Run algorithms — `BLRunner.py`

Runs one or more GRN inference algorithms on the specified datasets.

```bash
python BLRunner.py -c config-files/Curated/VSC.yaml
```

Each algorithm's output is written to `outputs/<dataset_dir>/<dataset_id>/<run_id>/<algorithm_id>/rankedEdges.csv`.

### 2. Evaluate results — `BLEvaluator.py`

Computes evaluation metrics by comparing each algorithm's ranked edge list to the ground truth network.

```bash
python BLEvaluator.py -c config-files/Curated/VSC.yaml [flags]
```

| Flag | Metric |
|------|--------|
| `-a` / `--auc` | AUPRC and AUROC |
| `-e` / `--epr` | Early precision ratio |
| `-s` / `--sepr` | Signed early precision (activation / inhibition) |
| `-r` / `--spearman` | Spearman correlation of predicted edge ranks |
| `-j` / `--jaccard` | Jaccard index of top-k predicted edges |
| `-t` / `--time` | Algorithm runtime |
| `-m` / `--motifs` | Network motif counts in top-k predicted networks |
| `-p` / `--paths` | Path length statistics on top-k predicted networks |
| `-b` / `--borda` | Borda-count edge aggregation across algorithms |

### 3. Plot results — `BLPlotter.py`

Generates publication-style figures from evaluation output.

```bash
python BLPlotter.py -c config-files/Curated/VSC.yaml -o ./plots [flags]
```

| Flag | Output file | Description |
|------|-------------|-------------|
| `-a` / `--auprc` | `AUPRC.pdf` | Precision-recall curves or box plots per dataset |
| `-r` / `--auroc` | `AUROC.pdf` | ROC curves or box plots per dataset |
| `-e` / `--epr` | `EPR.pdf` | Box plot of early precision values per algorithm |
| `--summary` | `Summary.pdf` | Heatmap of median AUPRC ratio and Spearman stability |
| `--epr-summary` | `EPRSummary.pdf` | Heatmap of AUPRC ratio, EPR ratio, and signed EPR ratios |
| `--all` | all of the above | Run all plots |

---

## Configuration

Config files are YAML and follow this structure:

```yaml
input_settings:
    input_dir: "inputs"
    dataset_dir: "Curated"
    datasets:
        - dataset_id: "mHSC-E"
          groundTruthNetwork: "GroundTruthNetwork.csv"
          runs:
            - run_id: "mHSC-E-500-1"
            - run_id: "mHSC-E-500-2"

    algorithms:
        - algorithm_id: "GENIE3"
          should_run: True
          params: {}

        - algorithm_id: "PPCOR"
          should_run: True
          params:
              pVal: 0.01

output_settings:
    output_dir: "outputs"
```

### `input_settings`

| Field | Required | Description |
|-------|----------|-------------|
| `input_dir` | Yes | Base directory containing all input datasets. Can be absolute or relative to the working directory. |
| `dataset_dir` | No | Subdirectory of `input_dir` that groups datasets by collection (e.g., `"Curated"`, `"Synthetic"`). |
| `datasets` | Yes | List of dataset groups. See **Dataset fields** below. |
| `algorithms` | Yes | List of algorithms to run. See **Algorithm fields** below. |

#### Dataset fields

| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `dataset_id` | Yes | — | Name of the dataset group. Used as a subdirectory under `dataset_dir`. |
| `should_run` | No | `[True]` | Set to `[False]` to skip this dataset entirely. |
| `groundTruthNetwork` | No | `GroundTruthNetwork.csv` | Filename of the ground truth edge list CSV, located in the dataset group directory (shared across all runs). |
| `runs` | Yes | — | List of individual run variants. See **Run fields** below. |

#### Run fields

Each entry under `runs` represents one replicate or condition variant. Input files are expected at `input_dir/dataset_dir/dataset_id/run_id/`.

| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `run_id` | Yes | — | Identifier for this run. Used as the subdirectory name within the dataset group directory. |
| `exprData` | No | `ExpressionData.csv` | Expression data filename, located in the run directory. |
| `pseudoTimeData` | No | `PseudoTime.csv` | Pseudotime data filename, located in the run directory. |

#### Algorithm fields

| Field | Required | Description |
|-------|----------|-------------|
| `algorithm_id` | Yes | Algorithm name. Must match one of the supported identifiers (see [Supported Algorithms](#supported-algorithms)). |
| `should_run` | Yes | Set to `True` to run this algorithm, `False` to skip it. |
| `params` | No | Dict of algorithm-specific parameters. Values are typically wrapped in a single-element list (e.g., `pVal: [0.01]`); the runner unwraps them automatically. |

### `output_settings`

| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `output_dir` | Yes | — | Base directory for all output files. Can be absolute or relative to the working directory. |
| `run_id` | No | — | When set, inserts an extra path segment between `output_dir` and `dataset_dir`. Useful for keeping outputs from separate experiment runs (e.g., different parameter sweeps) in the same base directory without overwriting each other. |

Output files are written to:
```
output_dir/[run_id/]dataset_dir/dataset_id/run_id/algorithm_id/rankedEdges.csv
```

---

## Supported Algorithms

| Algorithm | Runner file |
|-----------|-------------|
| GENIE3 | `BLRun/genie3Runner.py` |
| GRISLI | `BLRun/grisliRunner.py` |
| GRNBOOST2 | `BLRun/grnboost2Runner.py` |
| GRNVBEM | `BLRun/grnvbemRunner.py` |
| JUMP3 | `BLRun/jump3Runner.py` |
| LEAP | `BLRun/leapRunner.py` |
| PIDC | `BLRun/pidcRunner.py` |
| PPCOR | `BLRun/ppcorRunner.py` |
| SCODE | `BLRun/scodeRunner.py` |
| SCRIBE | `BLRun/scribeRunner.py` |
| SCSGL | `BLRun/scsglRunner.py` |
| SINCERITIES | `BLRun/sinceritiesRunner.py` |
| SINGE | `BLRun/singeRunner.py` |

---

## Citation

If you use BEELINE in your research, please cite:

Pratapa, A., Jalihal, A.P., Law, J.N., Bharadwaj, A., Murali, T.M. (2020) "Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data." _Nature Methods_, 17, 147–154.

- Publication: [https://www.nature.com/articles/s41592-019-0690-6](https://www.nature.com/articles/s41592-019-0690-6)
- Preprint: [https://doi.org/10.1101/642926](https://doi.org/10.1101/642926)

## Related Resources

- Input datasets: [https://doi.org/10.5281/zenodo.3378975](https://doi.org/10.5281/zenodo.3378975)
- BoolODE (synthetic data generator): [https://github.com/Murali-group/BoolODE](https://github.com/Murali-group/BoolODE)
- Docker images: [https://hub.docker.com/u/grnbeeline](https://hub.docker.com/u/grnbeeline)

---

## Directory Structure

```
.
├── BLRunner.py          # Entry point: run algorithms
├── BLEvaluator.py       # Entry point: evaluate results
├── BLPlotter.py         # Entry point: generate plots
├── initialize.sh        # Pull or build Docker images
├── setupAnacondaVENV.sh # Create/update BEELINE conda environment
├── BLRun/               # Algorithm runner classes
├── BLEval/              # Evaluation metric implementations
├── BLPlot/              # Plot generation implementations
├── config-files/        # YAML configuration files
├── inputs/              # Input datasets
└── outputs/             # Algorithm outputs (mirrors inputs/ structure)
```
