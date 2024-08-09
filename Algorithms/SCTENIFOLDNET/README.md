*This README.md file was generated on 2/20/2023 by Yiqi Su (yiqisu@vt.edu)*

**We would like to acknowledge professor Daniel Osorio for sharing the code for SCTENIFOLDNET.**

<!-- remove all comments (like this) before final save  -->

# scTenifoldNet: A Machine Learning Workflow for Constructing and Comparing Transcriptome-wide Gene Regulatory Networks from Single-Cell Data

This is the instruction on how to integrate the new GRN method SCTENIFOLDNET ([[Paper](https://doi.org/10.1016/j.patter.2020.100139)] [[GitHub](https://github.com/jamesjcai/ScTenifoldNet.jl)]) to BEELINE. 
Please follow the following steps:

1. **Create SCTENIFOLDNET folder:** Create a folder called SCTENIFOLDNET under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.

2. **Create runSCTENIFOLDNET.py script:** In the SCTENIFOLDNET folder, create an R script runSCTENIFOLDNET.r to learn graphs from target datasets. 

3. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. 

        FROM r-base:4.0.2

        LABEL maintainer = "Daniel Osorio <dcosorioh@tamu.edu>"

        USER root

        WORKDIR /

        RUN R -e "install.packages('https://cran.r-project.org/src/contrib/remotes_2.4.2.tar.gz', type = 'source')"

        RUN R -e "remotes::install_cran(pkgs = 'scTenifoldNet', quiet = TRUE)"

        RUN R -e "remotes::install_cran(pkgs = 'reshape2', quiet = TRUE)"

        RUN R -e "library(scTenifoldNet)"

        RUN R -e "library(reshape2)"

        COPY runSCTENIFOLDNET.R /

        RUN mkdir data/

        RUN apt-get update && apt-get install -y time


The Dockerfile will run the script runSCTENIFOLDNET.py within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for sctenifoldnet.

        cd $BASEDIR/Algorithms/SCTENIFOLDNET/
        docker build -q -t sctenifoldnet:base .
        if ([[ "$(docker images -q sctenifoldnet:base 2> /dev/null)" != "" ]]); then
            echo "Docker container for SCTENIFOLDNET is built and tagged as sctenifoldnet:base"
        else
            echo "Oops! Unable to build Docker container for SCTENIFOLDNET"
        fi

5. **Create sctenifoldnetRunner.py script:** After buliding the Docker image, create a Python script called sctenifoldnetRunner.py in Beeline/BLRun folder to setup a BLRun object so that it is able to read inputs and run sctenifoldnet inside the Docker image, and also parse the output for evaluation. Specifically, the sctenifoldnetRunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the input data file (i.e., expression data), and processes it into the format required by SCTENIFOLDNET. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of the input data file (i.e., expression data). It also specifies where the outputs are written. The docker container runs SCTENIFOLDNET when the parameters are passed. 
   - ``parseOutput()`` : This function reads the SCTENIFOLDNET-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

6. **Add SCTENIFOLDNET to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to SCTENIFOLDNET. 

    - add "import BLRun.sctenifoldnetRunner as SCTENIFOLDNET"
    - add "'SCTENIFOLDNET':SCTENIFOLDNET.generateInputs" to InputMapper
    - add "'SCTENIFOLDNET':SCTENIFOLDNET.run" to AlgorithmMapper
    - add "'SCTENIFOLDNET':SCTENIFOLDNET.parseOutput" to OutputParser

7. **Add SCTENIFOLDNET to config.yaml:** The final step is to add the new algorithm SCTENIFOLDNET and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.

        - name: "SCTENIFOLDNET"
          params:
              should_run: [True]
