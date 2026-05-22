# BEELINE Tutorial — Local Execution in GitHub Codespaces

This tutorial walks through running BEELINE against the BoolODE synthetic
networks **without Docker**, using four algorithms that are fast enough to run
directly in conda environments: **PEARSON**, **PPCOR**, **LEAP**, and
**GRNBOOST2**.

The goal is a self-contained end-to-end run that fits inside a free GitHub
Codespace.

---

## 1. Open the repo in Codespaces

From the GitHub repo page: **Code → Codespaces → Create codespace on
`ismb-tutorial`**. The default 4-core / 16 GB machine is sufficient.

The devcontainer at `.devcontainer/devcontainer.json` installs miniconda and
runs `utils/setupAnacondaVENV.sh` automatically the first time the codespace
boots, so step 2 below is already done by the time you attach. Expect the
initial create to take roughly 10-15 minutes for the four conda environments
to solve.

Alternatively, clone locally and run the setup yourself:

```bash
git clone -b ismb-tutorial https://github.com/Murali-group/Beeline.git
cd Beeline
```

---

## 2. Create the conda environments

A single shell script creates four conda environments — one shared BEELINE
environment plus one minimal env per algorithm that needs a non-Python runtime:

```bash
bash utils/setupAnacondaVENV.sh
```

| Environment | Purpose | Approx. size |
|---|---|---|
| `BEELINE` | Top-level pipeline (BLRunner / BLEvaluator / BLPlotter). Also runs PEARSON. | ~1 GB |
| `BEELINE-PPCOR` | R 4.4 + `r-ppcor` from conda-forge | ~350 MB |
| `BEELINE-LEAP` | R 4.2 + `r-leap` | ~350 MB |
| `BEELINE-GRNBOOST2` | Python 3.8 + `arboreto` from bioconda | ~500 MB |

The per-algorithm environments exist because PPCOR and LEAP need different R
versions, and `arboreto` is only on bioconda. The runners call them via
`conda run -n BEELINE-<ALG> ...` — no `conda activate` required.

Activate the main env so you can invoke the BL* scripts:

```bash
conda activate BEELINE
```

> **TODO:** Verify each env solves cleanly on a fresh Codespace. Pin versions
> in the three `utils/environment-*.yml` files if any drift.

---

## 3. Synthetic dataset

The 54 dataset directories the tutorial config exercises (6 BoolODE networks
x 3 cell counts x 3 replicates, ~3 MB total) are committed directly under
`inputs/Synthetic/`. There is nothing to download.

```
inputs/Synthetic/dyn-BF/dyn-BF-100/dyn-BF-100-1/{ExpressionData,GroundTruthNetwork,PseudoTime}.csv
inputs/Synthetic/dyn-BF/dyn-BF-100/dyn-BF-100-2/...
... and so on for dyn-BFC, dyn-CY, dyn-LI, dyn-LL, dyn-TF
```

The full v4 dataset (10 replicates per cell count, all curated scRNA-seq
networks) lives on Zenodo at
[doi.org/10.5281/zenodo.3378975](https://doi.org/10.5281/zenodo.3378975) if
you want to extend the tutorial beyond what's in this branch.

---

## 4. Run the algorithms

The tutorial config wires up all four local algorithms against every
synthetic network (6 networks × 3 cell counts × 3 replicates = 54 datasets):

```bash
python BLRunner.py -c config-files/tutorial.yaml
```

Output lands at `outputs/Tutorial/<network>/<run_id>/<algorithm_id>/rankedEdges.csv`.

> **TODO:** Time a full run end-to-end on a 4-core Codespace. If it exceeds
> ~30 minutes, prune `tutorial.yaml` further (e.g. drop the 5000-cell variants,
> or keep only `dyn-LI` and `dyn-BF`).

---

## 5. Evaluate the results

Compute AUPRC, AUROC, and Early Precision Ratio against the ground-truth
networks:

```bash
python BLEvaluator.py -c config-files/tutorial.yaml -a -e
```

Output: `outputs/Tutorial/<network>/-AUPRC.csv`, `-AUROC.csv`, `-EPR.csv`.

---

## 6. Plot the results

Generate publication-style figures:

```bash
python BLPlotter.py -c config-files/tutorial.yaml -o plots --all
```

You should see:

- `plots/AUPRC/<network>-AUPRC.{pdf,png}` — box plots across the 3 replicates
- `plots/AUROC/<network>-AUROC.{pdf,png}`
- `plots/EPR/<network>-EPR.{pdf,png}`
- `plots/Summary.{pdf,png}` — heatmap of median AUPRC / Spearman across all
  networks and algorithms

> **TODO:** Decide whether the tutorial should walk through interpreting the
> Summary heatmap (e.g. point out that PEARSON often outperforms more complex
> methods on the synthetic data — a key finding of the Pratapa et al. paper).

---

## What's next

- **Try a curated dataset.** Real scRNA-seq networks from
  `config-files/Curated/` (mHSC, mDC, hHEP, hESC) are harder problems. They
  require additional algorithms that don't run cleanly without Docker (SCODE,
  SINCERITIES, PIDC, etc.).
- **Add your own algorithm.** Subclass `BLRun.runner.Runner`, implement
  `generateInputs`, `run`, and `parseOutput`, then register the runner in
  `BLRun/runner_factory.py` (or wherever the dispatch lives).

---

## Open items (tasking record)

Use this list to track what still needs writing or verification before this
tutorial is ready to publish:

- [x] `.devcontainer/devcontainer.json` for one-click Codespaces boot
- [x] Confirm and document Zenodo v4 download URL + extraction layout
      (replaced by committing the tutorial subset directly under `inputs/`)
- [ ] End-to-end timing on a 4-core Codespace; prune dataset list if slow
- [x] Smoke-test `BEELINE-PPCOR` env builds and `runPPCOR.R` runs via `conda run`
- [x] Smoke-test `BEELINE-LEAP` env builds and `runLeap.R` runs via `conda run`
- [x] Smoke-test `BEELINE-GRNBOOST2` env builds and `runArboreto.py` runs via `conda run`
- [x] Verify `BLEvaluator.py -a -e` works against the `image: local` outputs
- [ ] Verify `BLPlotter.py --all` produces the expected files for this config
- [ ] Decide whether to include an interpretation section for the Summary plot
- [ ] (Stretch) Cache the synthetic dataset inside the devcontainer image so
      step 3 becomes a no-op
