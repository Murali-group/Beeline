*This README.md file was generated on 1/11/2023 by Yiqi Su (yiqisu@vt.edu)*
<!-- remove all comments (like this) before final save  -->

# GRNVBEM: A Bayesian framework for the inference of gene regulatory networks from time and pseudo-time series data

This is the instruction on how to integrate GRNVBEM ([[Paper](https://doi.org/10.1093/bioinformatics/btx605)] [[GitHub](https://github.com/mscastillo/GRNVBEM)]) to BEELINE. Please follow the following steps:

1. **Create GRNVBEM folder:** Create a folder called GRNVBEM under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.

2. **Create VBEM folder:** In the GRNVBEM folder, create a VBEM folder to store necessary files to run GRNVBEM:

    - ``readme.txt`` : instructions on GRNVBEM executation
    - ``run_GRNVBEM.sh`` : shell script for temporarily setting environment variables and executing the application
    - ``requiredMCRProducts.txt`` : required product number of MATLAB Compiler Runtime (MCR):
        - ``35000`` : MATLAB Runtime - Core         
        - ``35010`` : MATLAB Runtime - Numerics                     
        - ``35119`` : MATLAB Runtime - Statistics and Machine Learning Toolbox Addin

    - ``mccExcludedFiles.log`` : the list of excluded files
    - ``VBEM`` : example datasets
    - ``tOut.txt`` : parent-child weight probability score

3. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. 

        FROM ubuntu:18.04
        
        LABEL Maintainer = "Aditya Pratapa <adyprat@vt.edu>"
        
        USER root
        
        RUN apt-get -qq update && apt-get -qq install -y unzip xorg wget curl libstdc++6
        
        RUN mkdir /mcr-install && \
            mkdir /opt/mcr && \
            cd /mcr-install && \
            wget -q http://ssd.mathworks.com/supportfiles/downloads/R2019a/Release/0/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2019a_glnxa64.zip && \
            cd /mcr-install && \
            unzip MATLAB_Runtime_R2019a_glnxa64.zip && \
            ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
            cd / && \
            rm -rf mcr-install
            
        RUN mkdir VBEM/
        
        COPY VBEM/ /VBEM/
        
        WORKDIR VBEM/
        
        ENV LD_LIBRARY_PATH /opt/mcr/v96/runtime/glnxa64:/opt/mcr/v96/bin/glnxa64
        
        RUN mkdir data/
        
        RUN apt-get install time


The Dockerfile will run the GRNVBEM algorithm within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for GRNVBEM.

        cd $BASEDIR/Algorithms/GRNVBEM/
        docker build -q -t grnvbem:base .
        if ([[ "$(docker images -q grnvbem:base 2> /dev/null)" != "" ]]); then
            echo "Do
            cker container for GRNVBEM is built and tagged as grnvbem:base"
        else
            echo "Oops! Unable to build Docker container for GRNVBEM"
        fi

5. **Create grnvbemRunner.py script:** After buliding the Docker image, create a Python script called grnvbemRunner.py in Beeline/BLRun folder to setup BLRun objects to read inputs and run GRNVBEM inside the Docker image, and also parse the output for evaluation. Specifically, the grnvbemRunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the two input data files (i.e., expression data and the pseudotime), and processes them into the format required by the GRNVBEM. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data files (i.e., expression data and the pseudotime). It also specifies where the outputs are written. The docker container run GRNVBEM when the parameters are passed. 
   - ``parseOutput()`` : This function reads the grnvbem-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

6. **Add GRNVBEM to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to GRNVBEM. 

    - add "import BLRun.grnvbemRunner as GRNVBEM"
    - add "'GRNVBEM':GRNVBEM.generateInputs" to InputMapper
    - add "'GRNVBEM':GRNVBEM.run" to AlgorithmMapper
    - add "'GRNVBEM':GRNVBEM.parseOutput" to OutputParser

7. **Add GRNVBEM to config.yaml:** The final step is to add the new algorithm GRNVBEM and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.

        - name: "GRNVBEM"
          params: 
              should_run: [True]