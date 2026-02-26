*This README.md file was generated on 1/19/2023 by Yiqi Su (yiqisu@vt.edu)*
<!-- remove all comments (like this) before final save  -->

# SCNS: a graphical tool for reconstructing executable regulatory networks from single-cell genomic data

This is the instruction on how to integrate SCNS ([[Paper](https://doi.org/10.1186/s12918-018-0581-y)] [[GitHub](https://github.com/swoodhouse/SCNS-GUI)]) to BEELINE. Please follow the following steps:

1. **Create SCNS folder:** Create a folder called SCNS under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.

2. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. 

        FROM meval/modeleval:latest
        
        LABEL maintainer="Aditya Pratapa <adyprat@vt.edu>"
        
        USER root
        ENV TZ=America/New_York
        ENV DEBIAN_FRONTEND=noninteractive 
        RUN apt-get update
        RUN apt-get install -y build-essential patch
        RUN apt-get install -y fsharp \
	            mono-xbuild
        RUN git clone https://github.com/Z3Prover/z3
        WORKDIR z3
        RUN git checkout d6df51951f4cdc95f0dfd3b1297d04a465d8f2ca 
        
        RUN python2 scripts/mk_make.py --python
        RUN cd build && make
        RUN cd build && make install
        WORKDIR ..
        
        RUN export PATH="$PATH:/z3/build/"
        
        RUN  git clone https://github.com/swoodhouse/SCNS-Toolkit
        WORKDIR SCNS-Toolkit
        RUN git checkout 27cb7a349f239d450a9571b270abc38b053ad6b2
        
        WORKDIR SynthesisEngine
        
        RUN xbuild SynthesisEngine.sln
        
        RUN cp bin/Release/* .
        
        #RUN mkdir outDir
        #RUN mono SynthesisEngine.exe cmpStates.csv cmpEdges.csv cmpParameters.csv cmp_initial_states.txt cmp_target_states.txt outDir
        #WORKDIR outDir/
        #RUN cat Cebpa.txt
        
        RUN apt-get install time


The Dockerfile will run the SCNS algorithm within the Docker container.

3. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for SCNS.

        cd $BASEDIR/Algorithms/SCNS/s
        docker build -q -t scns:base .
        if ([[ "$(docker images -q scns:base 2> /dev/null)" != "" ]]); then
            echo "Do
            cker container for SCNS is built and tagged as scns:base"
        else
            echo "Oops! Unable to build Docker container for SCNS"
        fi

4s. **Create scnsRunner.py script:** After buliding the Docker image, create a Python script called scnsRunner.py in Beeline/BLRun folder to setup BLRun objects to read inputs and run SCNS inside the Docker image, and also parse the output for evaluation. Specifically, the scnsRunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the input data files (i.e., expression data and pseudotime data), and processes them into the format required by the SCNS. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data files (i.e., expression data and pseudotime data). It also specifies where the outputs are written. The docker container run SCNS when the parameters are passed. 
   - ``parseOutput()`` : This function reads the scns-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

5. **Add SCNS to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to SCNS. 

    - add "import BLRun.scnsRunner as SCNS"
    - add "'SCNS':SCNS.generateInputs" to InputMapper
    - add "'SCNS':SCNS.run" to AlgorithmMapper
    - add "'SCNS':SCNS.parseOutput" to OutputParser

6. **Add SCNS to config.yaml:** The final step is to add the new algorithm SCNS and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.

        - name: "SCNS"
          params: 
              should_run: [True]

