.. _boolode:

BoolODE
=======

Overview
--------

The BoolODE package converts a Boolean model to an ODE model and carries out stochastic
simulations. It reads a model definition file, and outputs the simulated expression data
for a given model.

.. code:: bash

   python src/BoolODE.py --path=data/test_vars.txt --max-time=5 --num-cells 10

.. toctree::
   :maxdepth: 2
   :caption: BoolODE Guide

   boolode-tutorial
   boolode-processing

.. _boolodeoptions:

Command Line Options
--------------------

For a full list of available command line options, run:

.. code:: bash

   python src/BoolODE.py --help

API Reference
-------------

.. note:: BoolODE is an external package. The following pages document its
   module structure; docstring content will only appear if the ``src/``
   package is present on the Python path at build time.

.. toctree::
   :maxdepth: 4

   src
