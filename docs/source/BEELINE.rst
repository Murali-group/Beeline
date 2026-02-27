.. _beeline:


BEELINE
=======

Overview
--------

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



Getting Started
---------------

The BEELINE pipeline interfaces with the implementations of various
algorithms through Docker containers.  Please follow this tutorial on
how to `install docker
<https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04>`_
on Ubuntu 18.04.



.. tip:: Setup docker to run without sudo using the following shell command

   .. code:: sh

                   sudo usermod -aG docker $USER

   See more details `here <https://askubuntu.com/questions/477551/how-can-i-use-docker-without-sudo>`_.

Once docker has been set up correctly, the next step is to create the
docker containers for each of the algorithms.  The script
``utils/initialize.sh`` pulls or builds these containers. Run the following in the
terminal

.. code:: bash

          . utils/initialize.sh

.. note:: This step will take a while!

We recommend using `Anaconda <https://www.anaconda.com>`_ for Python. Run the following command to automatically create an Anaconda virtual environment named BEELINE and install necessary libraries

.. code:: bash

    . utils/setupAnacondaVENV.sh


To compute proposed reconstructions using the 14 GRN algorithms on the example dataset, run

.. code:: python

          python BLRunner.py --config config-files/config.yaml

To compute areas under the ROC and PR curves for the proposed reconstructions, run

.. code:: bash

          python BLEvaluator.py --config config-files/config.yaml --auc

To display the complete list of evaluation options, run

.. code:: bash

          python BLEvaluator.py --help

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   beeline-tutorial
   algorithms
   beeline-developer
