In order to process the output of BoolODE simulations, configure settings under ``post_processsing``
in the config file. Ensure that ``do_post_processing`` is set to ``True``.

Generate Samples from Simulations
#################################

Carry out dimensionality reduction
##################################

These simulated trajectories can now In order to visualize the entire
dataset, we can carry out dimensionality reduction using t-SNE. A
script like the one below is a good starting point for this.

Here, the time point of each cell is inferred from the cell ID, and this information
is used to color each cell in the scatter plot. Darker colors imply early time points in the
simulation. The output should look like the following

.. figure:: figs/tree.png
   :align: center

   t-SNE visualization of BoolODE output
           
Generate dropouts from expression data
######################################

In order to mimic real single-cell expression datasets, BoolODE
includes ``genDropouts.py`` which implements dropouts as described in
the paper, by dropping expression values below ``DROP_CUTOFF`` using a
probability of ``DROP_PROB``. This script samples ``NCELLS`` from the
columns in the expression dataset ``EXPR``, and will throw an error if
the number is greater than the number of columns.



.. attention:: Ensure the ``--dropout`` option is passed. If not, ``genDropouts``
               will still randomly sample cells but not dropout any values.

Compute Pseudotime using Slingshot
##################################

.. note:: runSlingshot.py requires the SlingShot docker file! Please make sure
          docker has been setup and the container built.

          
Slingshot computes pseudotime trajectories for a given dataset by
first carrying out dimensionality reduction, then carrying out
*k*-means clustering on the low dimensional embedding in order to
compute trajectories. The number of clusters expected depend on the
features of the dataset. For instance, a dataset with two steady state
clusters should be specified with ``nClusters: 3``, specifying an
additional initial state cluster.


.. note:: ``runslingshot.py`` requires a 'pseudotime' file passed
           using the ``--pseudo`` option. This file should contain the
           actual simulation time from BoolODE, which can then be used
           to compare the quality of the inferred trajectory with the
           actual simulation time values. This file is NOT required by
           Slingshot itself.
               
