*This README.md file was generated on 12/10/2022 by Yiqi Su (yiqisu@vt.edu)*
<!-- remove all comments (like this) before final save  -->

# scSGL: Kernelized Signed Graph Learning for Single-Cell Gene Regulatory Network Inference

This is the instruction on how to integrate the new GRN method scSGL ([[Paper](https://doi.org/10.1093/bioinformatics/btac288)] [[GitHub](https://github.com/Single-Cell-Graph-Learning/scSGL)]) to BEELINE. Please follow the following steps:


1. **Create SCSGL folder:** Create a folder called SCSGL under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.

2. **Clone GitHub repositories:** Create a local copy of the scSGL GitHub repository using ` git clone https://github.com/Single-Cell-Graph-Learning/scSGL.git` and store it under the new method folder Beeline/Algorithms/SCSGL. 

3. **Create run_scSGL.py script:** In the SCSGL folder, create a Python script run_scSGL.py to learn signed graphs from target datasets. The arguments are as follows:

   - ``--expression_file`` : path of gene expression file
   - ``--ref_net_file`` : path of reference network file
   - ``--pos_density`` : the parameter controls the density of positive part of the learned signed graph with default as 0.45
   - ``--neg_density`` : the parameter controls the density of negative part of the learned signed graph with default as 0.45
   - ``--assoc`` : the association of signed graph with default as "correlation"
   - ``--out_file`` : path of output file

4. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. 



        FROM python:3.8  #This is the base image upon which necessary libraries are installed

        LABEL Maintainer="Yiqi Su <yiqisu@vt.edu>" #Additional information

        USER root #Set main user as root to avoid permission issues

        WORKDIR / #Sets current working directory

        COPY run_scSGL.py /  #Copy the main Python script which will perform necessary GRN computations

        COPY scSGL /scSGL  #Copy the original scSGL repo files stored in SCSGL folder 

        RUN apt-get update && apt-get install -y r-base time  #Installl time command to compute time taken

        UN pip install -r /scSGL/requirements.txt #Install the requirments and install R to conda environment

        RUN Rscript -e "install.packages('pcaPP')" #Install pcaPP to use zero inflated Kendall tau as a kernel

        RUN mkdir data/ #Make a directory to mount the folder containing input files

The Dockerfile will run the script run_scSGL.py within the Docker container.

5. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for scSGL.



        cd $BASEDIR/Algorithms/SCSGL/
        docker build -q -t scsgl:base .
        if ([[ "$(docker images -q scsgl:base 2> /dev/null)" != "" ]]); then
            echo "Docker container for SCSGL is built and tagged as scsgl:base"
        else
            echo "Oops! Unable to build Docker container for SCSGL"
        fi

6. **Create scsglRunner.py script:** After buliding the Docker image, create a Python script called scsglRunner.py in Beeline/BLRun folder to setup a BLRun object so that it is able to read inputs and run scSGL inside the Docker image, and also parse the output for evaluation. Specifically, the scsglRunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the two input data files (i.e., expression data and the ref_net file), and processes them into the format required by scSGL. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data files (i.e., expression data and the ref_net file), positive density (pos_density), negative density (neg_density), and association type (assoc). It also specifies where the outputs are written. The docker container runs scSGL when the parameters are passed. 
   - ``parseOutput()`` : This function reads the scSGL-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

7. **Add scSGL to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to scSGL. 

    - add "import BLRun.scsglRunner as SCSGL"
    - add "'SCSGL':SCSGL.generateInputs" to InputMapper
    - add "'SCSGL':SCSGL.run" to AlgorithmMapper
    - add "'SCSGL':SCSGL.parseOutput" to OutputParser

In addition, add "self.trueEdges = params['trueEdges']" in def __init__(self, params) within class Runner(object) for evaluation.


8. **Add scSGL to config.yaml:** The final step is to add the new algorithm scSGLC and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.



        - name: "SCSGL"
          params:
              should_run: [True]
              pos_density: [0.45]
              neg_density: [0.45]
              assoc: ["correlation"]
