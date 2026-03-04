Reproducing BEELINE Results
===========================

This page describes how to reproduce the results from the original BEELINE paper
(Pratapa et al., *Nature Methods*, 2020).


Environment Setup
-----------------

Before running any experiments, complete the standard BEELINE setup described in
:doc:`getting-started`:

1. Install and configure Docker.
2. Create the BEELINE conda environment:

   .. code:: bash

      . utils/setupAnacondaVENV.sh

3. Pull the algorithm Docker images:

   .. code:: bash

      . utils/initialize.sh


Downloading Data from Zenodo
-----------------------------

The datasets used in the paper are available on Zenodo at
`https://doi.org/10.5281/zenodo.3378975 <https://doi.org/10.5281/zenodo.3378975>`_.

.. note::

   **BEELINE v1.1 requires Zenodo dataset version v4.** Earlier versions of
   BEELINE should use Zenodo dataset version v3. If you are using an older
   version of BEELINE, please refer to the documentation for that version rather
   than this page. The significant change is a streamlined dircetory structure 
   that allows users to easily replicate certain figures.

The archive contains two top-level items:

- ``BEELINE-data/`` — the expression, pseudotime, and ground truth files for the
  Synthetic and Curated datasets, as well as the raw scRNA-seq experimental data.
- ``BEELINE-Networks/`` — reference networks (ChIP-seq and STRING) and transcription
  factor lists used for preprocessing the experimental datasets.

**Extracting the data**

Copy the contents of ``BEELINE-data/inputs/`` into the ``inputs/`` directory at the
root of this repository, and place ``BEELINE-Networks/`` at the repository root as
well:

.. code:: text

   inputs/
   ├── Curated/
   │   ├── GSD/            (GroundTruthNetwork.csv + 10 run subdirectories)
   │   ├── GSD-q50/        (50% dropout variant)
   │   ├── GSD-q70/        (70% dropout variant)
   │   ├── HSC/
   │   ├── HSC-q50/
   │   ├── HSC-q70/
   │   ├── mCAD/
   │   ├── mCAD-q50/
   │   ├── mCAD-q70/
   │   ├── VSC/
   │   ├── VSC-q50/
   │   └── VSC-q70/
   ├── Synthetic/
   │   ├── dyn-BF/
   │   │   ├── dyn-BF-100/ (10 run subdirectories)
   │   │   ├── dyn-BF-200/
   │   │   ├── dyn-BF-500/
   │   │   ├── dyn-BF-2000/
   │   │   └── dyn-BF-5000/
   │   ├── dyn-BFC/        (same cell-count structure)
   │   ├── dyn-CY/
   │   ├── dyn-LI/
   │   ├── dyn-LL/
   │   └── dyn-TF/
   └── scRNA-Seq/
       ├── hESC/
       ├── hHep/
       ├── mDC/
       ├── mESC/
       ├── mHSC-E/
       ├── mHSC-GM/
       └── mHSC-L/

   BEELINE-Networks/
   ├── Networks/
   │   ├── human/      (ChIP-seq and STRING networks)
   │   └── mouse/      (ChIP-seq and STRING networks)
   ├── human-tfs.csv
   └── mouse-tfs.csv


Synthetic and Curated Datasets (Figures 2 and 4)
-------------------------------------------------

Figures 2 and 4 evaluate algorithm performance on BoolODE-simulated (Synthetic)
and curated biological (Curated) datasets.

Each of the four Curated datasets (GSD, HSC, mCAD, VSC) is split into three
top-level directories by dropout level: full cell count (e.g. ``GSD``), 50%
dropout (``GSD-q50``), and 70% dropout (``GSD-q70``). Each directory contains a
shared ``GroundTruthNetwork.csv`` and 10 run subdirectories (e.g. ``GSD-2000-1``
through ``GSD-2000-10``), each with ``ExpressionData.csv`` and ``PseudoTime.csv``.

Each of the six Synthetic networks (dyn-BF, dyn-BFC, dyn-CY, dyn-LI, dyn-LL,
dyn-TF) contains five cell-count subdirectories (100, 200, 500, 2000, 5000
cells), each with 10 replicate run directories (e.g.
``dyn-BF/dyn-BF-100/dyn-BF-100-1/``). Each leaf run directory contains
``ExpressionData.csv``, ``PseudoTime.csv``, and ``GroundTruthNetwork.csv``.

**Step 1: Run the inference algorithms**

.. note::

   SCNS is not included in BEELINE v1.1 due to its long run time.

Config files for each dataset are in ``config-files/Curated/`` and
``config-files/Synthetic/``. Run ``BLRunner.py`` once for each config file:

.. code:: bash

   python BLRunner.py --config config-files/Curated/GSD.yaml
   python BLRunner.py --config config-files/Curated/HSC.yaml
   python BLRunner.py --config config-files/Curated/mCAD.yaml
   python BLRunner.py --config config-files/Curated/VSC.yaml
   python BLRunner.py --config config-files/Synthetic/dyn-BF.yaml
   python BLRunner.py --config config-files/Synthetic/dyn-BFC.yaml
   python BLRunner.py --config config-files/Synthetic/dyn-CY.yaml
   python BLRunner.py --config config-files/Synthetic/dyn-LI.yaml
   python BLRunner.py --config config-files/Synthetic/dyn-LL.yaml
   python BLRunner.py --config config-files/Synthetic/dyn-TF.yaml

**Step 2: Evaluate algorithm output**

Use ``BLEvaluator.py`` with the plot config files, which aggregate results across
all datasets for a given collection:

.. code:: bash

   # Synthetic networks
   python BLEvaluator.py --config config-files/Synthetic/PlotCuratedNetworks.yaml \
       --auc --spearman --epr --sepr

   # Curated networks
   python BLEvaluator.py --config config-files/Curated/PlotSimulatedNetworks.yaml \
       --auc --spearman --epr --sepr

**Step 3: Generate plots**

*Figure 2* — AUPRC ratio and Spearman stability summary heatmap:

.. code:: bash

   python BLPlotter.py --config config-files/Synthetic/PlotCuratedNetworks.yaml \
       --output outputs/plots/Synthetic --summary

   python BLPlotter.py --config config-files/Curated/PlotSimulatedNetworks.yaml \
       --output outputs/plots/Curated --summary

*Figure 4* — EPR and signed EPR summary heatmap:

.. code:: bash

   python BLPlotter.py --config config-files/Synthetic/PlotCuratedNetworks.yaml \
       --output outputs/plots/Synthetic --epr-summary

   python BLPlotter.py --config config-files/Curated/PlotSimulatedNetworks.yaml \
       --output outputs/plots/Curated --epr-summary

Output PDFs are written to the specified ``--output`` directories as
``Summary.pdf`` and ``EPRSummary.pdf`` respectively.


Experimental Datasets (Figures 5 and 6)
-----------------------------------------

Figures 5 and 6 include evaluation on seven real experimental scRNA-seq datasets
from ``inputs/scRNA-Seq/``:

.. csv-table::
   :widths: 20, 80

   "**hESC**", "Human embryonic stem cells"
   "**hHep**", "Human mature hepatocytes"
   "**mDC**", "Mouse dendritic cells"
   "**mESC**", "Mouse embryonic stem cells"
   "**mHSC-E**", "Mouse hematopoietic stem cells — erythroid lineage"
   "**mHSC-GM**", "Mouse hematopoietic stem cells — granulocyte/monocyte lineage"
   "**mHSC-L**", "Mouse hematopoietic stem cells — lymphoid lineage"

Each dataset directory contains ``ExpressionData.csv``, ``PseudoTime.csv``, and
``GeneOrdering.csv`` (a ranked list of genes with p-values and variance statistics).
These datasets do **not** include a pre-built ground truth network file; reference
networks must be sourced from ``BEELINE-Networks/``.

**Preprocessing with generateExpInputs.py**

The script ``utils/generateExpInputs.py`` produces a filtered ``ExpressionData.csv``
by selecting the most variable and statistically significant genes. It can
optionally restrict a reference network to only those genes. Run it once per
dataset/configuration you want to evaluate. For example, for hESC with 500
non-TF genes against the ChIP-seq reference network:

.. code:: bash

   python utils/generateExpInputs.py \
       --expFile inputs/scRNA-Seq/hESC/ExpressionData.csv \
       --geneOrderingFile inputs/scRNA-Seq/hESC/GeneOrdering.csv \
       --TFFile BEELINE-Networks/human-tfs.csv \
       --netFile BEELINE-Networks/Networks/human/hESC-ChIP-seq-network.csv \
       --pVal 0.01 \
       --numGenes 500 \
       --outPrefix inputs/Experimental/hESC/hESC-500

For mouse datasets, substitute ``BEELINE-Networks/mouse-tfs.csv`` and the
appropriate file from ``BEELINE-Networks/Networks/mouse/``. Run
``python utils/generateExpInputs.py --help`` for the full list of options.

The output files should be placed in the appropriate subdirectory of ``inputs/``
and referenced from a BEELINE config file as described in :ref:`configfiles`.

**Note on reproducibility**

Figures 5 and 6 are the result of a large number of runs across all seven
experimental datasets with significantly customized configurations. Fully
replicating these figures requires several weeks of compute time, and step-by-step
instructions for reproducing them exactly are not available. A user familiar with
BEELINE can configure the pipeline to reproduce these results from the
data available on Zenodo.
