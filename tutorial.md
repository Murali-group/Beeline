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

## 3. Verify the install

Run the verification script before kicking off any algorithms. It exercises
each conda env and confirms the tutorial data is in place:

```bash
bash utils/verify-install.sh
```

Expected output: 11 `PASS` lines and a `0 failed` summary. The script checks:

- `conda` is on PATH and all four envs exist
- `BEELINE` can import `pandas`, `pyyaml`, `networkx`, `scikit-learn`
- `BLRunner.py --help` runs (catches BL* import errors)
- `BEELINE-PPCOR` loads `library(ppcor)`
- `BEELINE-LEAP` loads `library(LEAP)`
- `BEELINE-GRNBOOST2` can `from arboreto.algo import grnboost2`
- The 6 synthetic networks each have a `GroundTruthNetwork.csv`
- `config-files/tutorial.yaml` parses

If anything fails, the script lists which check failed so you know where to
look. Don't move on until all checks pass.

---

## 4. Synthetic dataset

The 42 dataset directories the tutorial config exercises (6 BoolODE networks
x 7 replicates at 500 cells, ~5 MB total) are committed directly under
`inputs/Synthetic/`. There is nothing to download.

```
inputs/Synthetic/dyn-BF/dyn-BF-500/dyn-BF-500-1/{ExpressionData,PseudoTime}.csv
inputs/Synthetic/dyn-BF/dyn-BF-500/dyn-BF-500-2/...
... and so on for dyn-BFC, dyn-CY, dyn-LI, dyn-LL, dyn-TF
```

`GroundTruthNetwork.csv` lives one level up at the dataset-group level
(e.g. `inputs/Synthetic/dyn-BF/dyn-BF-500/GroundTruthNetwork.csv`) and is
shared across all replicates within that group.

The full v4 dataset (10 replicates per cell count, all curated scRNA-seq
networks) lives on Zenodo at
[doi.org/10.5281/zenodo.3378975](https://doi.org/10.5281/zenodo.3378975) if
you want to extend the tutorial beyond what's in this branch.

---

## 5. Run the algorithms

The tutorial config wires up all four local algorithms against every
synthetic network (6 networks × 7 replicates at 500 cells = 42 runs):

```bash
python BLRunner.py -c config-files/tutorial.yaml
```

Output lands at `outputs/Tutorial/<network>/<run_id>/<algorithm_id>/rankedEdges.csv`.

---

## 6. Evaluate the results

Compute AUPRC, AUROC, Early Precision Ratio, Spearman stability, and signed
EPR against the ground-truth networks. The Spearman (`-r`) and signed-EPR
(`--sepr`) outputs feed the summary heatmaps in step 7 — running with just
`-a -e` works but the heatmaps will be missing those rows.

```bash
python BLEvaluator.py -c config-files/tutorial.yaml -a -e -r --sepr
```

Output: per-dataset-group CSVs at
`outputs/Tutorial/<network>/<cell_count>/{AUPRC,AUROC,EarlyPrecision,EarlyPrecisionActivation,EarlyPrecisionInhibitory,Spearman}.csv`.

---

## 7. Plot the results

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

## 8. Try it on real scRNA-seq data

The synthetic pipeline above is a controlled setting: the ground-truth
network *is* what generated the expression. Real scRNA-seq is messier —
technical dropout, biological noise, and a ground truth that is only
partially known. This step runs the same four algorithms against a real
dataset so you can see how the numbers shift.

**Dataset:** 1071 Erythroid-lineage cells from Nestorowa et al. 2016
([GEO GSE81682](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81682)),
restricted to the 11 transcription factors of the Krumsiek 2011 HSC
Boolean circuit (Cebpa, Eklf, EgrNab, Fli1, Fog1, Gata1, Gata2, Gfi1, Pu1,
Scl, cJun). Expression values are `log2(CPM + 1)` computed directly from
GSE81682 raw HTSeq counts. Pseudotime is Nestorowa's slingshot output
along the Erythroid curve.

**Ground truth:** the 30-edge Krumsiek Boolean circuit, shipped with
BEELINE at `data/BEELINE-data/inputs/Curated/HSC/`.

**About EgrNab:** the Krumsiek circuit models `EgrNab` as a fused Boolean
node standing in for Egr1, Egr2, and Nab2. Since there is no `EgrNab` gene
to measure, this dataset uses **Egr1** alone as the representative — the
choice is data-driven: Egr1 has by far the highest across-cell variance
(6.75) versus Egr2 (0.27), Nab2 (4.39), or their mean composite (1.25).

Committed at `inputs/scRNA-Seq/mHSC-HSC/`:

```
inputs/scRNA-Seq/mHSC-HSC/
    GroundTruthNetwork.csv         # 30 edges from Krumsiek 2011
    Erythroid/
        ExpressionData.csv         # 11 genes x 1071 cells, log2(CPM+1)
        PseudoTime.csv             # 1071 cells, Erythroid pseudotime
```

Run the pipeline against a separate config that points at this dataset:

```bash
python BLRunner.py    -c config-files/tutorial-scrna.yaml
python BLEvaluator.py -c config-files/tutorial-scrna.yaml -a -e -r --sepr
python BLPlotter.py   -c config-files/tutorial-scrna.yaml -o plots --all
```

Output lands under `outputs/Tutorial-scRNA/mHSC-HSC/`. Expect the AUPRC
and EPR scores to be materially lower than on the synthetic networks —
that gap is the story. The Krumsiek circuit is a very simplified model of
Erythroid differentiation, and the ChIP-seq / Boolean edges it encodes
may or may not be reflected as co-expression relationships in the
Nestorowa data.

Two things to expect in the output:

- `Spearman.csv` will be empty for this dataset. Spearman stability is a
  cross-replicate metric, and mHSC-HSC has a single run (Erythroid) so
  there is nothing to compare against. This is not an error.
- GRNBOOST2's `output.txt` typically contains a `distributed.comm.core.CommClosedError`
  traceback near the end. This is dask cleanup noise; the ranked edge
  file is written correctly before it fires.

> **Note:** the tutorial ships the derived dataset directly; you do not
> need to re-download GSE81682 or re-derive the expression matrix. If you
> want to see how the file was built, look at `.work/nestorowa/` (not
> committed) or the derivation notes in this section.

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
- [x] Verify `BLPlotter.py --all` produces the expected files for this config
- [ ] Decide whether to include an interpretation section for the Summary plot
- [ ] (Stretch) Cache the synthetic dataset inside the devcontainer image so
      step 4 becomes a no-op
- [x] Add real-scRNA-seq step (Nestorowa 2016 mHSC + Krumsiek curated HSC network)
      at `inputs/scRNA-Seq/mHSC-HSC/` with `config-files/tutorial-scrna.yaml`
- [x] Smoke-test the four local algorithms against `tutorial-scrna.yaml`
      (PEARSON/PPCOR/LEAP/GRNBOOST2 all produce 110-121 ranked edges;
      BLEvaluator emits AUPRC/AUROC/EPR/signed-EPR cleanly; Spearman.csv
      is empty because mHSC-HSC has one run, not a bug)
- [ ] Verify `BLPlotter.py --all` works against the mHSC-HSC output
      (heatmaps with a single dataset column may render awkwardly)
- [ ] Consider whether the `EgrNab -> Egr1` substitution should be documented
      inline in tutorial.md's step 8 or in a separate methods appendix
