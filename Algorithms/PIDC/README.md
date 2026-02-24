*This README.md file was generated on 1/17/2023 by Yiqi Su (yiqisu@vt.edu)*
<!-- remove all comments (like this) before final save  -->

# PIDC: Inferring Gene Regulatory Network (GRN) using partial information decomposition and context

This is the instruction on how to integrate PIDC ([[Paper](https://doi.org/10.1016/j.cels.2017.08.014)] [[GitHub](https://github.com/hmutpw/PIDC)]) to BEELINE. Please follow the following steps:

1. **Create PIDC folder:** Create a folder called PIDC under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.

2. **Create installPackages.jl and runPIDC.jl scripts:** In the PIDC folder, create Julia scipts installPackages.jl and runPIDC.jl to install required packages and learn graph from target datasets, respectively.

3. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. 

        FROM julia:1.1.0-stretch
        
        LABEL maintainer="Aditya Pratapa <adyprat@vt.edu>"
        
        USER root
        
        WORKDIR /
        
        COPY installPackages.jl /
        
        # Julia libs we want
        
        RUN julia installPackages.jl
        
        COPY runPIDC.jl /
        
        RUN apt-get update && apt-get install time

The Dockerfile will run the PIDC algorithm within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for PIDC.

        cd $BASEDIR/Algorithms/PIDC/
        docker build -q -t pidc:base .
        if ([[ "$(docker images -q pidc:base 2> /dev/null)" != "" ]]); then
            echo "Do
            cker container for PIDC is built and tagged as pidc:base"
        else
            echo "Oops! Unable to build Docker container for PIDC"
        fi

5. **Create pidcRunner.py script:** After buliding the Docker image, create a Python script called pidcRunner.py in Beeline/BLRun folder to setup BLRun objects to read inputs and run PIDC inside the Docker image, and also parse the output for evaluation. Specifically, the pidcRunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the input data file (i.e., expression data), and processes them into the format required by the PIDC. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data file (i.e., expression data). It also specifies where the outputs are written. The docker container run PIDC when the parameters are passed. 
   - ``parseOutput()`` : This function reads the pidc-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

6. **Add PIDC to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to PIDC. 

    - add "import BLRun.pidcRunner as PIDC"
    - add "'PIDC':PIDC.generateInputs" to InputMapper
    - add "'PIDC':PIDC.run" to AlgorithmMapper
    - add "'PIDC':PIDC.parseOutput" to OutputParser

7. **Add PIDC to config.yaml:** The final step is to add the new algorithm PIDC and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.

        - name: "PIDC"
          params: 
              should_run: [True]