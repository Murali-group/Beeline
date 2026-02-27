.. _algorithms:

Supported Algorithms
====================

The following table lists the algorithms and the parameters they accept,
along with default parameter values.

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - **Algorithm**
     - **Input Parameters**
   * - SINCERITIES
     - - ``nBins`` : (Default = 10)
   * - SCODE
     - - ``z`` : (Default = 10)
       - ``nIter`` : (Default = 1000)
       - ``nRep`` : (Default = 6)
   * - SINGE
     - - ``lambda`` : (Default = 0.01)
       - ``dT`` : (Default = 15)
       - ``num_lags`` : (Default = 5)
       - ``kernel_width`` : (Default = 0.5)
       - ``prob_zero_removal`` : (Default = 0)
       - ``prob_remove_samples`` : (Default = 0.0)
       - ``family`` : (Default = 'gaussian')
       - ``num_replicates`` : (Default = 6)
   * - PPCOR
     - - ``pVal`` : p-value cutoff (Default = 0.01)
   * - PIDC
     - None
   * - PEARSON
     - None
   * - LEAP
     - - ``maxLag`` : (Default = 0.33)
   * - JUMP3
     - None
   * - SCRIBE
     - - ``delay`` : (Default = 5)
       - ``method`` : Any of 'RDI', 'uRDI', 'cRDI', or 'ucRDI' (Default = 'ucRDI')
       - ``lowerDetectionLimit`` : (Default = 0)
       - ``expressionFamily`` : If mRNA counts, use 'negbinomial' (Default = 'uninormal')
       - ``log`` : Log transform expression values (Default = False)
       - ``ignorePT`` : Ignore pseudotime (Default = True)
   * - SCSGL
     - - ``pos_density`` : Density of positive edges in the inferred network
       - ``neg_density`` : Density of negative edges in the inferred network
       - ``assoc`` : Association measure used for graph learning
   * - GRNVBEM
     - None
   * - GRISLI
     - - ``L`` : (Default = 10)
       - ``R`` : (Default = 3000)
       - ``alphaMin`` : (Default = 0.0)
   * - GENIE3
     - None
   * - GRNBOOST2
     - None
