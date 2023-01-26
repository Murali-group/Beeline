*This README.md file was generated on 1/10/2023 by Yiqi Su (yiqisu@vt.edu)*
<!-- remove all comments (like this) before final save  -->

# Arboreto

Arboreto (https://arboreto.readthedocs.io/en/latest/index.html) currently supports two gene regulatory network (GRN) inference algorithms: GRNBoost2 ([[Paper](https://doi.org/10.1093/bioinformatics/bty916)] [[GitHub](https://github.com/aertslab/GRNBoost)]) and GENIE3 ([[Paper](https://doi.org/10.1371/journal.pone.0012776)] [[GitHub](https://github.com/vahuynh/GENIE3)]).

This is the instruction on how to integrate Arboreto to BEELINE. Please follow the following steps:

1. **Create ARBORETO folder:** Create a folder called ARBORETO under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.


2. **Create runArboreto.py script:** In the ARBORETO folder, create a Python script runArboreto.py to learn graphs from target datasets. The arguments are as follows:

   - ``--algo`` : specify either GRNBoost2 or GENIE3 to run
   - ``--inFile`` : path of input tab-separated expression file
   - ``--outFile`` : file where the output network is stored 


3. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. 


        FROM continuumio/anaconda3:2018.12
        
        LABEL Maintainer="Aditya Pratapa <adyprat@vt.edu>"
        
        USER root
        
        RUN apt-get update
        
        RUN conda install -y -c bioconda/label/cf201901 arboreto=0.1.5 pandas=0.24.0
        
        COPY runArboreto.py /
        
        RUN mkdir data/
        
        RUN apt-get install time

The Dockerfile will run the script runArboreto.py within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for Arboreto.


        cd $BASEDIR/Algorithms/ARBORETO/
        docker build -q -t arboreto:base .
        if ([[ "$(docker images -q arboreto:base 2> /dev/null)" != "" ]]); then
            echo "Docker container for ARBORETO is built and tagged as arboreto:base"
        else
            echo "Oops! Unable to build Docker container for ARBORETO"
        fi

5. **Create grnboost2Runner.py and genie3Runner.py scripts:** After buliding the Docker image, create Python scripts called grnboost2Runner.py and genie3Runner.py in Beeline/BLRun folder to setup BLRun objects to read inputs and run GRNBoost2 and GENIE3 respectively inside the Docker image, and also parse the output for evaluation. Specifically, both the two Python scripts contain three functions:

   - ``generateInputs()`` : This function reads the expression data file, and processes it into the format required by the corresponding algorithm. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data file. It also specifies where the outputs are written. The docker container run the corresponding algorithm when the parameters are passed. 
   - ``parseOutput()`` : This function reads the algorithm-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

6. **Add GRNBoost2 and GENIE3 to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to GRNBoost2 and GENIE3, respectively. 

    For GRNBoost2:
    - add "import BLRun.grnboost2Runner as GRNBOOST2"
    - add "'GRNBOOST2':GRNBOOST2.generateInputs" to InputMapper
    - add "'GRNBOOST2':GRNBOOST2.run" to AlgorithmMapper
    - add "'GRNBOOST2':GRNBOOST2.parseOutput" to OutputParser

    For GENIE3:
    - add "import BLRun.genie3Runner as GENIE3"
    - add "'GENIE3':GENIE3.generateInputs" to InputMapper
    - add "'GENIE3':GENIE3.run" to AlgorithmMapper
    - add "'GENIE3':GENIE3.parseOutput" to OutputParser


7. **Add GRNBoost2 and GENIE3 to config.yaml:** The final step is to add the new algorithms GRNBoost2 and GENIE3 and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.


        For GRNBoost2:
        - name: "GRNBOOST2"
          params: 
              should_run: [True]

        For GENIE3:
        - name: "GENIE3"
          params: 
              should_run: [True]


