*This README.md file was generated on 1/19/2023 by Yiqi Su (yiqisu@vt.edu)*
<!-- remove all comments (like this) before final save  -->

# SINGE: ANetwork inference with Granger causality ensembles on single-cell transcriptomics

This is the instruction on how to integrate SINGE ([[Paper](https://doi.org/10.1016/j.celrep.2022.110333)] [[GitHub](https://github.com/gitter-lab/SINGE)]) to BEELINE. Please follow the following steps:

1. **Create SINGE folder:** Create a folder called SINGE under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.

2. **Download files:** In the SINGE folder, download sh files including SINGe.sh, run_SINGE_Aggregate.sh, and run_SINGE_GLG_Test.sh from GitHub (https://github.com/gitter-lab/SINGE). In addition, download tests folder and keep the enviroment.yaml and run_SINGE_Test.sh files. 

3. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. 

        # MATLAB R2018a in a Debian environment
        FROM amarburg/matlab-runtime:R2018a
        
        USER root
        
        RUN apt-get update && \
            apt-get -y install libxt6 bzip2 time octave && \
            rm -rf /var/lib/apt/lists/*
            
        # Install Miniconda3 following https://hub.docker.com/r/continuumio/miniconda3/dockerfile
        # Python is only needed for testing SINGE and could be removed from the base
        # to reduce the image size
        ENV PATH /opt/conda/bin:$PATH
        RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O ~/miniconda.sh && \
            /bin/bash ~/miniconda.sh -b -p /opt/conda && \
            rm ~/miniconda.sh && \
            /opt/conda/bin/conda clean --all && \
            ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
            echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
            echo "conda activate base" >> ~/.bashrc && \
            conda init bash
        
        ENV SINGE_ROOT /usr/local/SINGE
        WORKDIR /usr/local/SINGE
        ENTRYPOINT ["/usr/local/SINGE/SINGE.sh", "/usr/local/MATLAB/MATLAB_Runtime/v94"]
        
        # Install conda test environment
        COPY tests/environment.yml tests/
        RUN conda env create -f tests/environment.yml

        # Copy wrapper scripts for compiled MATLAB executables
        COPY SINGE.sh .
        COPY run_SINGE_GLG_Test.sh .
        COPY run_SINGE_Aggregate.sh .
        COPY tests/run_SINGE_Test.sh tests/

        ENV SINGE_RUNNING_IN_DOCKER 1

        # Download the compiled SINGE executables from the stable release
        # md5sum of v0.4.1 SINGE_GLG_Test is c50ec7bc13e287eca340c9d19d8bc27d/
        RUN tag=v0.4.1 && \
            wget --quiet https://github.com/gitter-lab/SINGE/releases/download/$tag/SINGE_Test && \
            wget --quiet https://github.com/gitter-lab/SINGE/releases/download/$tag/SINGE_GLG_Test && \
            wget --quiet https://github.com/gitter-lab/SINGE/releases/download/$tag/SINGE_Aggregate && \
            wget --quiet https://github.com/gitter-lab/SINGE/releases/download/$tag/code.md5 && \
            chmod u+x SINGE_* && \
            mv SINGE_Test tests/SINGE_Test

        # Download an intermediate version of the compiled SINGE executables for testing
        # Download the md5sums of the source .m files and binaries
        #RUN md5=c50ec7bc13e287eca340c9d19d8bc27d && \
        #    wget --quiet https://www.biostat.wisc.edu/~gitter/tmp/$md5/SINGE_Test && \
        #    wget --quiet https://www.biostat.wisc.edu/~gitter/tmp/$md5/SINGE_GLG_Test && \
        #    wget --quiet https://www.biostat.wisc.edu/~gitter/tmp/$md5/SINGE_Aggregate && \
        #    wget --quiet https://www.biostat.wisc.edu/~gitter/tmp/$md5/code.md5 && \
        #    chmod u+x SINGE_* && \
        #    mv SINGE_Test tests/SINGE_Test


The Dockerfile will run the SINGE algorithm within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for SINGE.

        cd $BASEDIR/Algorithms/SINGE/
        docker build -q -t singe:base .
        if ([[ "$(docker images -q singe:base 2> /dev/null)" != "" ]]); then
            echo "Do
            cker container for SINGE is built and tagged as singe:base"
        else
            echo "Oops! Unable to build Docker container for SINGE"
        fi

5. **Create singeRunner.py script:** After buliding the Docker image, create a Python script called singeRunner.py in Beeline/BLRun folder to setup BLRun objects to read inputs and run SINGE inside the Docker image, and also parse the output for evaluation. Specifically, the singeRunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the input data files (i.e., expression data and pseudotime data), and processes them into the format required by the SINGE. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data files (i.e., expression data and pseudotime data). It also specifies where the outputs are written. The docker container run SINGE when the parameters are passed. 
   - ``parseOutput()`` : This function reads the singe-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

6. **Add SINGE to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to SINGE. 

    - add "import BLRun.singeRunner as SINGE"
    - add "'SINGE':SINGE.generateInputs" to InputMapper
    - add "'SINGE':SINGE.run" to AlgorithmMapper
    - add "'SINGE':SINGE.parseOutput" to OutputParser

7. **Add SINGE to config.yaml:** The final step is to add the new algorithm SINGE and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.

        - name: "SINGE"
          params: 
              should_run: [True]
              lambda: [0.01]
              dT: [15]
              num_lags: [5]
              kernel_width: [0.5]
              prob_zero_removal: [0]
              prob_remove_samples: [0.0]
              family: ["gaussian"]
              num_replicates: [6]