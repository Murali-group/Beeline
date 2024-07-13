# Implementation of the different methods for GRN inference.

The scripts provided are:

- `genie3.R`: Functions to infer gene regulatory networks with GENIE3 with or without refinement via chromatin accessibility data.
- `l0l2.R`: Functions to infer gene regulatory networks with L0L2 regression with or without refinement via chromatin accessibility data.
- `mutual_info.R`: Functions to infer gene regulatory networks with Mutual Information with (i.e. MICA) or without refinement via chromatin accessibility data.
- `spearman.R`: Functions to infer gene regulatory networks with Spearman correlation with or without refinement via chromatin accessibility data.
- `rnd_predictor.R`: Functions to produce a random gene regulatory networks with or without refinement via chromatin accessibility data.


# MICA : Mutual Information Chromatin Accessibility

MICA attempts to use information from RNA seq and Chromatin Accessibility information to improve GRN prediction. 
It uses four methods for the same : Spearman Correlation, GENIE3, Mutual Information and L0L2 Sparse Regression


MICA ([[Paper](https://doi.org/10.1101/2023.02.03.527081)] [[GitHub](https://github.com/SydneyBioX/scTIE)]) 
This is the instruction on how to integrate MICA to BEELINE. Please follow the following steps:

1. **Create MICA folder:** Create a folder called MICA under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implementations.


2. **Create runMICA.R script:** In the MICA folder, create an R script runMICA.py to learn graphs from target datasets. The arguments are as follows:

   - ``--method`` : specify either spearman, l0l2, genie3, mica to run
   - ``--expressionFile`` : path of input gene expression file
   - ``--atacFile`` : path to the chromatin accessibility and peaks files
   - ``--regFile`` : path to file with information about regulators 


3. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom.

```
FROM r-base:4.3.1

RUN apt-get update && apt-get install -y

USER root

WORKDIR /

RUN mkdir data/

COPY runMICA.R /
COPY spearman.R /
COPY l0l2.R /
COPY genie3.R /
COPY mutual_info.R /



LABEL maintainer="Vishnu Madhav <vm.ibab@gmail.com>"
LABEL version="1.0"
LABEL description="Mutual Information Chromatin Accessibility : Gene Regulatory network building"
```
The Dockerfile will run the script runMICA.R within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for MICA.

cd $BASEDIR/Algorithms/MICA/
docker build -q -t mica:base .
if ([ $? = 0 ] && [[ "$(docker images -q mica:base 2> /dev/null)" != "" ]]); then
  echo "Docker container for MICA is built and tagged as mica:base"
elif [ "$(docker images -q mica:base 2> /dev/null)" != "" ]; then
    echo "Docker container failed to build, but an existing image exists at mica:base"
else
    echo "Oops! Unable to build Docker container for MICA"
fi

5. **Create micaRunner.py scripts:** After building the Docker image, create Python scripts called micaRunner.py in Beeline/BLRun folder to set up BLRun objects to read inputs and run MICA inside the Docker image, and also parse the output for evaluation. Specifically, they contain three functions:

   - ``generateInputs()`` : This function reads the expression data file, and processes it into the format required by the corresponding algorithm. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data file. It also specifies where the outputs are written. The docker container run the corresponding algorithm when the parameters are passed. 
   - ``parseOutput()`` : This function reads the algorithm-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

6. **Add GMICA to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to MICA, respectively. 

    For GRNBoost2:
    - add "import BLRun.micaRunner as MICA"
    - add "'MICA':MICA.generateInputs" to InputMapper
    - add "'MICA':MICA.run" to AlgorithmMapper
    - add "'MICA':MICA.parseOutput" to OutputParser


7. **Add MICA to config.yaml:** The final step is to add the new algorithms MICA and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.


        For MICA:
        - name: "MICA"
          params: 
              should_run: [True]
