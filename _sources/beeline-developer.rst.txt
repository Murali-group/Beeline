.. _blrunguide:

BLRun Details
#############

The ``BLRun.py`` creates a :obj:`BLRun` object for each algorithm
specified in the config file. Each :obj:`BLRun` object should contain three
modules

1. ``generateInputs()`` : This function reads the three input data
   files, and processes them into the format required by the given
   algorithm
2. ``run()`` : A function to construct a system command with the
   appropriate command line parameters to be passed to the docker
   container in order to run a given algorithm
3. ``parseOutput()`` : A function to read the algorithm-specific
   output and reformat it into a standard format

The evaluation scripts in the final step of the pipeline expect the
inferred networks from each algorithm to be a comma-separated file
with the following format:

.. code:: text

          Gene1,Gene2,EdgeWeight
          reg1,targ1,edgeweight

where the first line are the column names, and the subsequent lines
contain the edges predicted by the network. The Gene1 column should
contain regulators, the Gene2 column the targets, and EdgeWeight
column the absolute value of the weight predicted for edge (regulator,
target).

.. _blevalguide:

BLEval Details
##############

This is the list of options of currently implemented evaluation functions

Command line arguments:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
                        Configuration file containing list of datasets
                        algorithms and output specifications.
  -a, --auc             Compute median of areas under Precision-Recall and ROC
                        curves. Calls :mod:`BLEval.computeAUC`.
  -j, --jaccard         Compute median Jaccard index of predicted top-k
                        networks for each algorithm for a given set of
                        datasets generated from the same ground truth network. Calls :mod:`BLEval.computeJaccard`.
  -r, --spearman        Compute median Spearman Corr. of predicted edges for
                        each algorithm for a given set of datasets generated
                        from the same ground truth network.  Calls :mod:`BLEval.computeSpearman`.
  -t, --time            Analyze time taken by each algorithm for a. Calls :mod:`BLEval.parseTime`.
  -e, --epr             Compute median early precision. Calls :mod:`BLEval.computeEarlyPrec`.
  -s, --sepr            Analyze median (signed) early precision for activation
                        and inhibitory edges. :mod:`BLEval.computeSignedEPrec`.
  -m, --motifs          Compute network motifs in the predicted top-k
                        networks. Calls :mod:`BLEval.computeNetMotifs`.

