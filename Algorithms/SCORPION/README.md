*This README.md file was generated on 2/4/2023 by Yiqi Su (yiqisu@vt.edu)*

**We would like to acknowledge professors Daniel Osorio, S. Stephen Yi and Marieke L. Kuijjer for sharing the code for SCORPION.**

<!-- remove all comments (like this) before final save  -->

# SCORPION: Single-Cell Oriented Reconstruction of PANDA (https://sites.google.com/a/channing.harvard.edu/kimberlyglass/tools/panda) Individually Optimized Gene Regulatory Network

This is the instruction on how to integrate the new GRN method SCORPION ([[Paper](https://doi.org/10.1101/2023.01.20.524974)] [[GitHub](https://github.com/kuijjerlab/SCORPION)]) to BEELINE. 
Please follow the following steps:

1. **Create SCORPION folder:** Create a folder called SCORPION under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.

2. **Create runSCORPION.py script:** In the SCORPION folder, create an R script runSCORPION.r to learn graphs from target datasets. 

3. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. 


        FROM r-base:4.2.0

        LABEL maintainer = "Daniel Osorio <daniel.osoriohurtado@childrens.harvard.edu>"

        USER root

        WORKDIR /

        RUN R -e "install.packages('https://cran.r-project.org/src/contrib/remotes_2.4.2.tar.gz', type = 'source')"

        RUN R -e "install.packages('reshape2')"

        RUN R -e "remotes::install_github('kuijjerlab/SCORPION')"

        RUN R -e "library(SCORPION)"

        RUN R -e "library(reshape2)"

        COPY runSCORPION.R /

        RUN mkdir data/

        RUN apt-get update && apt-get install time

The Dockerfile will run the script runSCORPION.py within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for scorpion.


        cd $BASEDIR/Algorithms/SCORPION/
        docker build -q -t scorpion:base .
        if ([[ "$(docker images -q scorpion:base 2> /dev/null)" != "" ]]); then
            echo "Docker container for SCORPION is built and tagged as scorpion:base"
        else
            echo "Oops! Unable to build Docker container for SCORPION"
        fi

5. **Create scorpionRunner.py script:** After buliding the Docker image, create a Python script called scorpionRunner.py in Beeline/BLRun folder to setup a BLRun object so that it is able to read inputs and run scorpion inside the Docker image, and also parse the output for evaluation. Specifically, the scorpionRunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the input data file (i.e., expression data), and processes it into the format required by SCORPION. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of the input data file (i.e., expression data). It also specifies where the outputs are written. The docker container runs SCORPION when the parameters are passed. 
   - ``parseOutput()`` : This function reads the SCORPION-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

6. **Add SCORPION to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to SCORPION. 

    - add "import BLRun.scorpionRunner as SCORPION"
    - add "'SCORPION':SCORPION.generateInputs" to InputMapper
    - add "'SCORPION':SCORPION.run" to AlgorithmMapper
    - add "'SCORPION':SCORPION.parseOutput" to OutputParser

7. **Add SCORPION to config.yaml:** The final step is to add the new algorithm SCORPION and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.


        - name: "SCORPION"
          params:
              should_run: [True]
