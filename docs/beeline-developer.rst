.. _blrunguide:

Please follow the following steps to add a new GRN inferenmence method to BEELINE. 

Adding a new GRN inference method
#################################

In order to ensure easy set-up and portability, all the GRN algorithms in BEELINE are provided within their own separate Docker instances. This also avoids conflicting libraries/software versions that may arise from the GRN algorithm implmentations. More resources on Docker can be found `here <https://docs.docker.com/get-started/resources/>`_. In order to add a new GRN method to BEELINE, you will need the following components.

1. **Create a Dockerfile:** To create a Docker image for your GRN algorithm locally, you will need to start with a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. `BEELINE <https://github.com/Murali-group/Beeline/tree/master/Algorithms/>`_ currently provides Dockerfiles for GRN implementations in `Python <https://github.com/Murali-group/Beeline/tree/master/Algorithms/ARBORETO>`_, `R <https://github.com/Murali-group/Beeline/tree/master/Algorithms/PPCOR>`_, `MATLAB <https://github.com/Murali-group/Beeline/blob/master/Algorithms/GRISLI/Dockerfile>`_, and `Julia <https://github.com/Murali-group/Beeline/tree/master/Algorithms/PIDC>`_. These Dockerfiles can be used as a template for adding other GRN algorithms. Note that in order to run MATLAB-based code, we recommend first builing a standalone MATLAB application for your GRN algorithm  `MATLAB Runtime <https://www.mathworks.com/help/compiler/create-and-install-a-standalone-application-from-matlab-code.html>`_ and then setting up your Dockerfile to simply execute the pre-compiled binaries. For example, the `Dockerfile <https://github.com/Murali-group/Beeline/blob/master/Algorithms/PPCOR/Dockerfile>`_ using R which runs the script `runPPCOR.R <https://github.com/Murali-group/Beeline/blob/master/Algorithms/PPCOR/runPPCOR.R>`_ within the Docker container is as follows:

.. code:: Dockerfile

    FROM r-base:3.5.3 #This is the base image upon which necessary libraries are installed

    LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>" #Additional information

    RUN apt-get update && apt-get install time #Installl time command to compute time taken
    
    USER root #Set main user as root to avoid permission issues

    WORKDIR / #Sets current working directory

    RUN R -e "install.packages('https://cran.r-project.org/src/contrib/ppcor_1.1.tar.gz', type = 'source')" #Installs a specific version of PPCOR package

    COPY runPPCOR.R / #Copy the main R script which will perform necessary GRN computations

    RUN mkdir data/ #Make a directory to mount the folder containing input files

The Dockerfile must then be placed under the `Algorithms/<algorithm-name> <https://github.com/Murali-group/Beeline/tree/master/Algorithms/>`_ folder in BEELINE. Then the following lines must be added to the `initialize.sh <https://github.com/Murali-group/Beeline/blob/master/initialize.sh>`_ script.

.. code:: bash

    cd $BASEDIR/Algorithms/<algorithm-name>/
    docker build -q -t <algorithm-name>:base .
    echo "Docker container for <algorithm-name> is built and tagged as <algorithm-name>:base"

Replace <algorithm-name> with the name of the new algorithm that is being added to BEELINE.

2. **Create a <algorithm-name>Runner.py script:** Once the Docker image is built locally using the above Dockerfile, we will then need to setup a :obj:`BLRun` object which will read necessary inputs, runs GRN algorithm inside the Docker image, and finally parses output so that evaluation can be performed using :obj:`BLEval`.
Most of the GRN algorithms in BEELINE require a gene-by-cell matrix provided as input, along with a pesudotime ordering of cells, and any additional manually specified parameters. These details are specified as inputs to BEELINE using :ref:`configfiles` as mentioned earlier. In our current implementation, we provide a separate csv file for expression matrix and pseudotime files, whereas algorithm parameters are specified as command line arguments. However, this can be modified easily to specify even the parameters using a separate file by simply providing its path in the config file. Each <algorithm-name>Runner.py script should contain the following three functions:

   - ``generateInputs()`` : This function reads the two input data files (i.e., expression data and the pseudotime), and processes them into the format required by the given algorithm. For example, the algorithm may only require the cells in the expression matrix to be ordered in pseudotime, instead of exact pesudotime values as input. Moreover, if the algorithm requires that you run the GRN inference on each trajectory separately, we can this function to write the separate expression matrices contianing only cells  of a particular trajectory. 
   - ``run()`` : This function constructs a "docker run" system command with the appropriate command line parameters to be passed to the docker container in order to run a given algorithm. The docker run command also needs to mount the folder containing the input files using the "-v" flag and specify where the outputs are written. If the algorithm needs to be run separately on each trajectory, it can be simply called in a loop inside this function.
   - ``parseOutput()`` : This function reads the algorithm-specific outputs and formats it into a ranked edgelist comma-separated file in the following format which can be subsequently used by :obj:`BLEval`

.. code:: text

          Gene1,Gene2,EdgeWeight
          reg1,targ1,edgeweight

where the first line are the column names, and the subsequent lines contain the edges predicted by the network. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). Note that in cases where the algorithm requires that you run the GRN inference on each trajectory separately, a single output network is obtained by finding the maximum score for a given edge across the GRNs computed for each individual trajectory.  Please ensure that all of the above three functions in your script accept arguments of type :obj:`BLRun.runner.Runner`. The <algorithm-name>Runner.py script must then be placed under the `BLRun/ <https://github.com/Murali-group/Beeline/tree/master/BLRun/>`_  folder in BEELINE.  

3. **Add the new alorithm to runner.py:** The next step is to integrate the new algorithm within :obj:`BLRun` object. This can be achieved by adding the above three modules from the above step, i.e, ``generateInputs()``, ``run()``, and ``parseOutput()`` to `runner.py <https://github.com/Murali-group/Beeline/blob/master/BLRun/runner.py>`_.

4. **Add the new alorithm to config.yaml:** The final step is to add the new algorithm and any necessary parameters to the cofig.yaml. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.


.. code:: text

      algorithms:
          - name: <algorithm-name>
            params:
                # Any parameters that can be passed to this algorithm
                paramX: [0.1]
                paramY: ['Text']
                should_run: [True] # or False


.. _blevalguide:

Adding a new evaluation technique
#################################

BEELINE currently supports several evaluation techniques, namely, area under ROC and PR curves, eary precision, early signed precision, time and memory consumption, network motifs, jaccard similarities, and spearman coefficients. A new evaluation technique can be easily added to BEELINE using the following procedure.

1. Add the script containing the new evaluation technique to the `BLEval/ <https://github.com/Murali-group/Beeline/tree/master/BLEval/>`_ folder.

2. The next step is to integrate the new technique into the :obj:`BLEval` object. This can be achieved by adding it as a new module under :obj:`BLEval` in the `BLEval/__init__.py <https://github.com/Murali-group/Beeline/blob/master/BLEval/__init__.py>`_ script. Please ensure that your script can accept arguments of type :obj:`BLEval`.

3. The final step is to add a command line option to perform the evaluation to `BLEvaluator.py <https://github.com/Murali-group/Beeline/blob/master/BLEvaluator.py>`_.


