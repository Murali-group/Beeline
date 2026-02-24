*This README.md file was generated on 1/12/2023 by Yiqi Su (yiqisu@vt.edu)*
<!-- remove all comments (like this) before final save  -->

# GRISLI: Gene Regulation Inference for Single-cell with LInear differential equations and velocity inference

This is the instruction on how to integrate GRISLI2 ([[Paper](https://doi.org/10.1093/bioinformatics/btaa576)] [[GitHub](https://github.com/PCAubin/GRISLI)]) to BEELINE. Please follow the following steps:

1. **Create GRISLI folder:** Create a folder called GRISLI under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.


2. **Create runGRISLI folder:** In the GRISLI folder, create a runGRISLI folder to store necessary files to run GRISLI:

    - ``spams-matlab-v2.6`` : the SPAMS toolbox (http://spams-devel.gforge.inria.fr/) to solve the lasso problem
    - ``readme.txt`` : instructions on GRISLI executation
    - ``run_GRISLI.sh`` : shell script for temporarily setting environment variables and executing the application
    - ``requiredMCRProducts.txt`` : required product number of MATLAB Compiler Runtime (MCR):
        - ``35000``	: MATLAB Runtime - Core			
        - ``35010``	: MATLAB Runtime - Numerics						
        - ``35119``	: MATLAB Runtime - Statistics and Machine Learning Toolbox Addin

    - ``mccExcludedFiles.log`` : the list of excluded files
    - ``GRISLI`` : example datasets


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
            
        RUN mkdir runGRISLI
        
        COPY runGRISLI/ /runGRISLI/
        
        WORKDIR runGRISLI/
        
        ENV LD_LIBRARY_PATH /opt/mcr/v96/runtime/glnxa64:/opt/mcr/v96/bin/glnxa64
        
        RUN mkdir -p /root/.mcrCache9.6/GRISLI0/GRISLI/
        
        RUN cp -r /runGRISLI/spams-matlab-v2.6/ /root/.mcrCache9.6/GRISLI0/GRISLI/
        
        RUN mkdir data/
        
        RUN apt-get update
        
        RUN apt-get install -y libgomp1 --fix-missing
        
        RUN apt-get install time

The Dockerfile will run the GRISLI algorithm within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for GRISLI.


        cd $BASEDIR/Algorithms/GRISLI/
        docker build -q -t grisli:base .
        if ([[ "$(docker images -q grisli:base 2> /dev/null)" != "" ]]); then
            echo "Docker container for GRISLI is built and tagged as grisli:base"
        else
            echo "Oops! Unable to build Docker container for GRISLI"
        fi

5. **Create grisliRunner.py script:** After buliding the Docker image, create a Python script called grisliRunner.py in Beeline/BLRun folder to setup BLRun objects to read inputs and run GRISLI inside the Docker image, and also parse the output for evaluation. Specifically, the grisliRunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the two input data files (i.e., expression data and the pseudotime), and processes them into the format required by the GRISLI. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data files (i.e., expression data and the pseudotime). It also specifies where the outputs are written. The docker container run GRISLI when the parameters are passed. 
   - ``parseOutput()`` : This function reads the grisli-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

6. **Add GRISLI to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to GRISLI. 

    - add "import BLRun.grisliRunner as GRISLI"
    - add "'GRISLI':GRISLI.generateInputs" to InputMapper
    - add "'GRISLI':GRISLI.run" to AlgorithmMapper
    - add "'GRISLI':GRISLI.parseOutput" to OutputParser


7. **Add GRISLI to config.yaml:** The final step is to add the new algorithm GRISLI and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.


        - name: "GRISLI"
          params: 
              should_run: [True]
              L: [5]
              R: [1500]
              alphaMin: [0.0]