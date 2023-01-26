*This README.md file was generated on 1/17/2023 by Yiqi Su (yiqisu@vt.edu)*
<!-- remove all comments (like this) before final save  -->

# LEAP: constructing gene co-expression networks for single-cell RNA-sequencing data using pseudotime ordering

This is the instruction on how to integrate LEAP ([[Paper](https://doi.org/10.1093/bioinformatics/btw729)] [[GitHub](https://github.com/cran/LEAP)]) to BEELINE. Please follow the following steps:

1. **Create LEAP folder:** Create a folder called LEAP under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.

2. **Create runLeap.R script:** In the LEAP folder, create an R scipt runLeap.R to learn graph from target datasets.

3. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. 

        FROM r-base
        
        LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"
        
        USER root
        
        WORKDIR /
        
        RUN R -e "install.packages('https://cran.r-project.org/src/contrib/LEAP_0.2.tar.gz', type = 'source')"
        
        COPY runLeap.R /
        
        RUN mkdir data/
        
        RUN apt-get update && apt-get install -y time

The Dockerfile will run the LEAP algorithm within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for LEAP.

        cd $BASEDIR/Algorithms/LEAP/
        docker build -q -t leap:base .
        if ([[ "$(docker images -q leap:base 2> /dev/null)" != "" ]]); then
            echo "Do
            cker container for LEAP is built and tagged as leap:base"
        else
            echo "Oops! Unable to build Docker container for LEAP"
        fi

5. **Create leapRunner.py script:** After buliding the Docker image, create a Python script called leapRunner.py in Beeline/BLRun folder to setup BLRun objects to read inputs and run LEAP inside the Docker image, and also parse the output for evaluation. Specifically, the leapRunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the two input data files (i.e., expression data and the pseudotime), and processes them into the format required by the LEAP. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data files (i.e., expression data and the pseudotime). It also specifies where the outputs are written. The docker container run LEAP when the parameters are passed. 
   - ``parseOutput()`` : This function reads the leap-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

6. **Add LEAP to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to LEAP. 

    - add "import BLRun.leapRunner as LEAP"
    - add "'LEAP':LEAP.generateInputs" to InputMapper
    - add "'LEAP':LEAP.run" to AlgorithmMapper
    - add "'LEAP':LEAP.parseOutput" to OutputParser

7. **Add LEAP to config.yaml:** The final step is to add the new algorithm LEAP and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.

        - name: "LEAP"
          params: 
              should_run: [True]
              # Default maxLag value is 0.33
              maxLag: [0.33]