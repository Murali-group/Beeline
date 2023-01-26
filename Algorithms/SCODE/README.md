*This README.md file was generated on 1/19/2023 by Yiqi Su (yiqisu@vt.edu)*
<!-- remove all comments (like this) before final save  -->

# SCODE: an efficient regulatory network inference algorithm from single-cell RNA-Seq during differentiation

This is the instruction on how to integrate SCODE ([[Paper](https://doi.org/10.1093/bioinformatics/btx194)] [[GitHub](https://github.com/hmatsu1226/SCODE)]) to BEELINE. Please follow the following steps:

1. **Create SCODE folder:** Create a folder called SCODE under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.

2. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. 

        FROM r-base
        
        LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"
        
        USER root
        
        WORKDIR /
        
        
        RUN apt-get update 
        
        # RUN add-apt-repository ppa:ubuntu-toolchain-r/test
        
        RUN apt-get install -y  gcc
        
        RUN apt-get install -y git
        
        RUN apt-get install -y ruby
        
        RUN git clone https://github.com/hmatsu1226/SCODE
        
        WORKDIR SCODE/
        
        RUN git checkout a0512f8ec29aac188c9c27a8e89ddd2464e6d84d
        
        RUN R -e "install.packages('https://cran.r-project.org/src/contrib/MASS_7.3-51.3.tar.gz', type = 'source')"
        
        # RUN ruby run_R.rb data/exp_train.txt data/time_train.txt out 100 4 356 100 2
        
        RUN apt-get install -y time


The Dockerfile will run the SCODE algorithm within the Docker container.

3. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for SCODE.

        cd $BASEDIR/Algorithms/SCODE/
        docker build -q -t scode:base .
        if ([[ "$(docker images -q scode:base 2> /dev/null)" != "" ]]); then
            echo "Do
            cker container for SCODE is built and tagged as scode:base"
        else
            echo "Oops! Unable to build Docker container for SCODE"
        fi

4. **Create scodeRunner.py script:** After buliding the Docker image, create a Python script called scodeRunner.py in Beeline/BLRun folder to setup BLRun objects to read inputs and run SCODE inside the Docker image, and also parse the output for evaluation. Specifically, the scodeRunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the input data files (i.e., expression data and pseudotime data), and processes them into the format required by the SCODE. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data files (i.e., expression data and pseudotime data). It also specifies where the outputs are written. The docker container run SCODE when the parameters are passed. 
   - ``parseOutput()`` : This function reads the scode-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

5. **Add SCODE to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to SCODE. 

    - add "import BLRun.scodeRunner as SCODE"
    - add "'SCODE':SCODE.generateInputs" to InputMapper
    - add "'SCODE':SCODE.run" to AlgorithmMapper
    - add "'SCODE':SCODE.parseOutput" to OutputParser

6. **Add SCODE to config.yaml:** The final step is to add the new algorithm SCODE and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.

        - name: "SCODE"
          params: 
              should_run: [True]
              z: [10]
              nIter: [1000]
              nRep: [6]