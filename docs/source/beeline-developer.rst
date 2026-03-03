Developer Guide
===============

BEELINE is designed to be extensible. This section describes how to add new
GRN inference methods and new evaluation techniques to the pipeline.

.. _blrunguide:

Adding a new GRN inference method
----------------------------------

Please follow the following steps to add a new GRN inference method to BEELINE.

In order to ensure easy set-up and portability, all the GRN algorithms in BEELINE are provided within their own separate Docker instances. This also avoids conflicting libraries/software versions that may arise from the GRN algorithm implementations. More resources on Docker can be found `here <https://docs.docker.com/get-started/resources/>`_. In order to add a new GRN method to BEELINE, you will need the following components.

1. **Create a Dockerfile:** To create a Docker image for your GRN algorithm locally, you will need to start with a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. `BEELINE <https://github.com/Murali-group/Beeline/tree/master/Algorithms/>`_ currently provides Dockerfiles for GRN implementations in `Python <https://github.com/Murali-group/Beeline/tree/master/Algorithms/ARBORETO>`_, `R <https://github.com/Murali-group/Beeline/tree/master/Algorithms/PPCOR>`_, `MATLAB <https://github.com/Murali-group/Beeline/blob/master/Algorithms/GRISLI/Dockerfile>`_, and `Julia <https://github.com/Murali-group/Beeline/tree/master/Algorithms/PIDC>`_. These Dockerfiles can be used as a template for adding other GRN algorithms. Note that in order to run MATLAB-based code, we recommend first building a standalone MATLAB application for your GRN algorithm using `MATLAB Runtime <https://www.mathworks.com/help/compiler/create-and-install-a-standalone-application-from-matlab-code.html>`_ and then setting up your Dockerfile to simply execute the pre-compiled binaries. For example, the `Dockerfile <https://github.com/Murali-group/Beeline/blob/master/Algorithms/PPCOR/Dockerfile>`_ using R which runs the script `runPPCOR.R <https://github.com/Murali-group/Beeline/blob/master/Algorithms/PPCOR/runPPCOR.R>`_ within the Docker container is as follows:

.. code:: Dockerfile

    FROM r-base:4.2.2

    LABEL maintainer="Aditya Pratapa <adyprat@vt.edu>"

    USER root

    WORKDIR /

    RUN R -e "install.packages('ppcor', repos='https://cloud.r-project.org', dependencies=TRUE)"

    COPY runPPCOR.R /

    RUN mkdir data/

    RUN apt-get update -o Acquire::AllowInsecureRepositories=true \
        && apt-get install -y --allow-unauthenticated debian-archive-keyring \
        && apt-get update \
        && apt-get install -y time

The Dockerfile must then be placed under the ``Algorithms/<algorithm-name>/`` folder in BEELINE. Then the following lines must be added to the `utils/initialize.sh <https://github.com/Murali-group/Beeline/blob/master/initialize.sh>`_ script.

.. code:: bash

    cd $BASEDIR/Algorithms/<algorithm-name>/
    docker build -q -t <algorithm-name>:base .
    echo "Docker container for <algorithm-name> is built and tagged as <algorithm-name>:base"

Replace ``<algorithm-name>`` with the name of the new algorithm being added to BEELINE.

2. **Create a <algorithm-name>Runner.py script:** Once the Docker image is built locally using the above Dockerfile, we will then need to set up a Runner class that reads necessary inputs, runs the GRN algorithm inside the Docker image, and finally parses the output so that evaluation can be performed using :obj:`BLEval`.
Most of the GRN algorithms in BEELINE require a gene-by-cell matrix provided as input, along with a pseudotime ordering of cells, and any additional manually specified parameters. These details are specified as inputs to BEELINE using :ref:`configfiles` as mentioned earlier. Each ``<algorithm-name>Runner.py`` script must define a class that inherits from :class:`BLRun.runner.Runner`, which is an Abstract Base Class (ABC). Concrete subclasses must implement all three abstract methods listed below, or Python will raise a ``TypeError`` at instantiation time:

   - ``generateInputs(self)`` : This method reads the two input data files (i.e., expression data and the pseudotime), and processes them into the format required by the given algorithm. Input files are read from ``self.input_dir`` and processed files are written to ``self.working_dir``.
   - ``run(self)`` : This method constructs a ``docker run`` system command with the appropriate command line parameters to run the algorithm inside its container. The working directory ``self.working_dir`` is mounted as ``/usr/working_dir`` inside the container.
   - ``parseOutput(self)`` : This method reads the algorithm-specific outputs and formats them into a ranked edgelist written to ``self.output_dir/rankedEdges.csv`` in the following format:

.. code:: text

          Gene1,Gene2,EdgeWeight
          reg1,targ1,edgeweight

where the first line contains the column names, and the subsequent lines contain the edges predicted by the network. The ``Gene1`` column should contain regulators, ``Gene2`` the targets, and ``EdgeWeight`` the absolute value of the weight predicted for edge (regulator, target). Each method must accept only ``self`` as its argument, matching the abstract method signatures defined in :class:`BLRun.runner.Runner`. The ``<algorithm-name>Runner.py`` script must then be placed under the `BLRun/ <https://github.com/Murali-group/Beeline/tree/master/BLRun/>`_ folder in BEELINE.

3. **Register the new algorithm in BLRunner.py:** The next step is to integrate the new algorithm with the pipeline. Import the new Runner class at the top of `BLRunner.py <https://github.com/Murali-group/Beeline/blob/master/BLRunner.py>`_ and add it to the ``RUNNERS`` dictionary:

.. code:: python

    from BLRun.<algorithm-name>Runner import <AlgorithmName>Runner

    RUNNERS = {
        ...
        '<ALGORITHM-NAME>': <AlgorithmName>Runner,
    }

4. **Add the new algorithm to the config file:** The final step is to add the new algorithm and any necessary parameters to the config YAML. Each algorithm entry requires an ``algorithm_id``, a Docker ``image`` tag, a ``should_run`` flag, and an optional ``params`` block.

.. code:: text

      algorithms:
          - algorithm_id: <algorithm-name>
            image: "<docker-image>:base"
            should_run: [True]  # or [False]
            params:
                # Any parameters that can be passed to this algorithm
                paramX: [0.1]
                paramY: ['Text']


.. _blevalguide:

Adding a new evaluation technique
----------------------------------

BEELINE currently supports several evaluation techniques, namely, area under ROC and PR curves, early precision, early signed precision, time and memory consumption, network motifs, Jaccard similarities, Spearman coefficients, Borda aggregation, and path length statistics. A new evaluation technique can be easily added to BEELINE using the following procedure.

1. Add the script containing the new evaluation technique to the `BLEval/ <https://github.com/Murali-group/Beeline/tree/master/BLEval/>`_ folder.

2. The next step is to integrate the new technique into the :obj:`BLEval` object. This can be achieved by adding it as a new module under :obj:`BLEval` in the `BLEval/__init__.py <https://github.com/Murali-group/Beeline/blob/master/BLEval/__init__.py>`_ script. Please ensure that your script can accept arguments of type :obj:`BLEval`.

3. The final step is to add a command line option to perform the evaluation to `BLEvaluator.py <https://github.com/Murali-group/Beeline/blob/master/BLEvaluator.py>`_.
