# TENET: Transfer Entropy Network

This is the instruction on how to integrate TENET ([[Paper](https://doi.org/10.1093/nar/gkaa1014)] [[GitHub](https://github.com/neocaleb/TENET)]) to BEELINE. Please follow the following steps:

1. **Create TENET folder:** Create a folder called TENET under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implementations.

2. **Clone the TENET repository:** In the TENET folder, clone the [[TENET GitHub repository](https://github.com/neocaleb/TENET)]. 

    - **Warning:** There is a known issue with running TENET's multi-threaded implementation, as well as some indexing issues with TENET's single core version. BEELINE uses TENET's single core version, which requires the following manual fix inside ``TENET/TENETsinglecore``:

        * To match Beeline's method of generating input data, modify line 46 by adding a ``[1:]`` to the end. This allows for correct indexing of the pandas DataFrame object when printing results.

        * On lines 107 and 108, modify the indexing of the TEmatrix object by removing all `-1`'s. 
    
    - This should result in the following lines:

        * **Line 46**: ``gene_name=line.replace('\n','').split(',')[1:]``

        * **Line 107**: ``TEmatrix[int(temp[0])][int(temp[1])]=float(temp[2])``

        * **Line 108**: ``TEmatrix[int(temp[1])][int(temp[0])]=float(temp[3])``

3. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. 


        FROM python:3

        LABEL Maintainer="Tim Wilson <timwilson@vt.edu>"

        USER root

        WORKDIR /

        RUN apt-get update && apt-get install -y time && \ 
            apt-get install wget -y &&\
            # install openmpi
            wget https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.3.tar.bz2 && \
            tar xf openmpi-5.0.3.tar.bz2 && \
            cd openmpi-5.0.3 && \
            mkdir build && \
            cd build && \
            ../configure && \
            make -j 4 && \
            apt-get install default-jre -y && \
            apt-get install default-jdk -y
            

        RUN pip install jpype1 statsmodels numpy    

        COPY TENET /TENET

        RUN mkdir /data

        WORKDIR /TENET

The Dockerfile will run the TENET algorithm within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to ``initialize.sh`` to create the Docker image for TENET.


        cd $BASEDIR/Algorithms/TENET/
        docker build -q -t tenet:base .
        if ([ $? = 0 ] && [[ "$(docker images -q tenet:base 2> /dev/null)" != "" ]]); then
            echo "Docker container for TENET is built and tagged as tenet:base"
        elif [ "$(docker images -q tenet:base 2> /dev/null)" != "" ]; then
            echo "Docker container failed to build, but an existing image exists at tenet:base"
        else
            echo "Oops! Unable to build Docker container for TENET"
        fi

5. **Create tenetRunner.py script:** After building the Docker image, create a Python script called ``tenetRunner.py`` in Beeline/BLRun folder to setup BLRun objects to read inputs and run TENET inside the Docker image, and also parse the output for evaluation. Specifically, the ``tenetRunner.py`` script contains three functions:

    - ``generateInputs()``: This function reads the two input data files (i.e., expression data and the pseudotime), and processes the into the format required by TENET.

    - ``run()``: This function constructs a "docker run" system command with parameters, including the path of the input data files (i.e., the expression data and the pseudotime). It also specifies where the outputs are written. The docker container runs TENET when the parameters are passed.

    - ``parseOutput()``: This function reads the TENET-specific output (i.e., ``outFile.txt``) and formats it into a ranked edgelist comma-separated file (i.e., ``rankedEdges.csv``) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator, target). The ranked edgelist file will be subsequently used by the BLEval object.

6. **Add TENET to runner.py:** Next, update ``runner.py`` in Beeline/BLRun folder by adding information related to TENET.

    - add "``import BLRun.tenetRunner as TENET``"
    - add "``'TENET':TENET.generateInputs``" to InputMapper
    - add "``'TENET':TENET.run``" to AlgorithmMapper
    - add "``'TENET':TENET.parseOutput``" to OutputParser

7. **Add TENET to config.yaml:** The final step is to add the new algorithm TENET and any necessary parameters to the ``config.yaml`` located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time even though multiple parameters can be passed onto the single parameter object. 

        - name: "TENET"
          params:
                should_run: [True]
                historyLength: [1]