.. BEELINE documentation master file, created by
   sphinx-quickstart on Fri Feb 27 19:31:37 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:github_url:  https://github.com/murali-group/Beeline

.. _beeline:

Welcome to BEELINE's documentation!
=====================================

.. figure:: ../figs/overview-graphic.png
   :align: center

   Overview of BEELINE

This site serves as the documentation for the code released under the
**Benchmarking gEnE reguLatory network Inference from siNgle-cEll
transcriptomic data** (BEELINE) project.  The first part of this
project is the BEELINE pipeline which provides tools to
evaluate the performance of algorithms for the reconstruction of gene
regulatory networks (GRNs) from single-cell RNAseq data. The second
part of this project is :ref:`BoolODE`, a tool to automatically
convert a Boolean model to an ODE model, and subsequently carry out
stochastic simulations.

BEELINE provides a set of tools for evaluating methods that infer gene
regulatory networks (GRN) from single-cell gene expression
data. The BEELINE framework is divided into the following modules:

* :ref:`blrun` : contains BEELINE's Runner module, a Python
  wrapper for 14 GRN inference algorithms with options to add new
  methods.
* :ref:`bleval` : contains BEELINE's Evaluation module that
  provides easy-to-use tools for evaluating GRN reconstructions.
* :ref:`blplot` : contains BEELINE's plotting module for
  visualizing output from BLEval.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Home <self>
   getting-started
   running-beeline
   reproducing-results
   algorithms
   beeline-developer
   BoolODE
   modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
