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

Tutorial
--------

.. include:: boolode-tutorial.rst

.. _geninputs:

Generating inputs for BEELINE
-----------------------------

.. include:: boolode-processing.rst

.. _boolodeoptions:

Command Line Options
--------------------

.. include:: usage.rst

API Reference
-------------

.. include:: boolode-reference.rst


