*This README.md file was generated on 1/20/2023 by Yiqi Su (yiqisu@vt.edu)*
<!-- remove all comments (like this) before final save  -->

# SCRIBE: Towards inferring causal gene regulatory networks from single cell expression Measurements

This is the instruction on how to integrate SCRIBE ([[Paper](https://doi.org/10.1101/426981)] [[GitHub](https://github.com/cole-trapnell-lab/Scribe)]) to BEELINE. Please follow the following steps:

1. **Create SCRIBE folder:** Create a folder called SCRIBE under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.

2. **Create runScribe.R script:** In the SCRIBE folder, create an R script runScribe.py for to learn graphs from target datasets. The arguments are as follows:

    - ``--expressionFile`` : path of gene expression file
    - ``--cellFile`` : path to comma separated file containing data on cells
    - ``--geneFile"`` :  path to comma separated file containing data on genes
    - ``--newCellDataSet`` : path to .RDS file containing an object of type newCellDataSet
    - ``--lowerDetectionLimit`` : single float value to pass as an argument for newCellDataSet function
    - ``--expressionFamily`` : VGAM family function name to be used for expression response variables
    - ``--method`` : method name for Scribe. Can be any one of 'RDI', 'cRDI', 'uRDI', or 'ucRDI'
    - ``--delay`` : comma separated list of delay values for Scribe
    - ``--log`` : Log transform expression values
    - ``--outPrefix`` : path to write output files. Required
    - ``--outFile`` : outFile name to write the output ranked edges
    - ``--ignorePT`` : ignores pseudotime computed using monocle and uses experiment time
   

3. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. 

        FROM bioconductor/release_base2:R3.5.3_Bioc3.8
        
        LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"
        
        USER root
        
        WORKDIR /
        
        Run R -e "BiocManager::install('devtools',version=3.8)"
        
        RUN apt-get update && apt-get install -y libhdf5-dev \
        libxml2-dev \
        libudunits2-dev \
        imagemagick \
        zlib1g-dev \
        libfreetype6-dev
        
        RUN DEBIAN_FRONTEND=noninteractive apt-get -y install xorg \
        libx11-dev \
        libglu1-mesa-dev 
        
        RUN apt-get install -y r-cran-rgl
        
        
        RUN R -e "BiocManager::install('HiveR',ref = 3.8)"
        
        RUN R -e "BiocManager::install(c('lattice','Matrix','irlba'),ref = 3.8, update = FALSE)"
        
        RUN R -e "BiocManager::install(c('cluster','dplyr'), ref = 3.8, update = FALSE)"
        
        # Needs monocle 2.8
        
        RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('monocle')" 
        
        RUN git clone https://github.com/cole-trapnell-lab/RANNinf
        
        WORKDIR RANNinf/
        
        RUN git checkout 9d4a3d781c8b74a01ec2018ca7137b22a0e583b6
        
        WORKDIR /
        
        RUN R CMD build RANNinf
        
        RUN R -e "BiocManager::install('RcppArmadillo', ref = 3.8, update = FALSE)"
        
        RUN R -e "install.packages('RANNinf_2.5.0.99.tar.gz', repo = NULL, type ='source')"
        
        RUN git clone https://github.com/cole-trapnell-lab/Scribe
        
        WORKDIR Scribe/
        
        RUN git checkout 4ba98500764adbce4a59be508d94b279bbfcfb31
        
        RUN R -e "BiocManager::install(c('cowplot','lpSolveAPI'), ref = 3.8, update = FALSE)"
        
        RUN R -e "install.packages('Scribe_0.1.tar.gz', repo = NULL, type ='source')"
        
        WORKDIR /
        
        RUN R -e "BiocManager::install(c('optparse'), ref = 3.8, update = FALSE)"
        
        RUN mkdir data/
        
        COPY runScribe.R /
        
        RUN apt-get install -y time

The Dockerfile will run the SCRIBE algorithm within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for SCRIBE.

        cd $BASEDIR/Algorithms/SCRIBE/
        docker build -q -t scribe:base .
        if ([[ "$(docker images -q scribe:base 2> /dev/null)" != "" ]]); then
            echo "Do
            cker container for SCRIBE is built and tagged as scribe:base"
        else
            echo "Oops! Unable to build Docker container for SCRIBE"
        fi

5. **Create scribeRunner.py script:** After buliding the Docker image, create a Python script called scribeRunner.py in Beeline/BLRun folder to setup BLRun objects to read inputs and run SCRIBE inside the Docker image, and also parse the output for evaluation. Specifically, the scribeRunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the input data files (i.e., expression data and pseudotime data), and processes them into the format required by the SCRIBE. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data files (i.e., expression data and pseudotime data). It also specifies where the outputs are written. The docker container run SCRIBE when the parameters are passed. 
   - ``parseOutput()`` : This function reads the scribe-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

6. **Add SCRIBE to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to SCRIBE. 

    - add "import BLRun.scribeRunner as SCRIBE"
    - add "'SCRIBE':SCRIBE.generateInputs" to InputMapper
    - add "'SCRIBE':SCRIBE.run" to AlgorithmMapper
    - add "'SCRIBE':SCRIBE.parseOutput" to OutputParser

7. **Add SCRIBE to config.yaml:** The final step is to add the new algorithm SCRIBE and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.

        - name: "SCRIBE"
          params: 
              should_run: [True]
              ### required parameters
              # a list of delay values
              delay: ["5"]
              # any of 'RDI', 'uRDI', 'cRDI', or 'ucRDI'
              method: ['ucRDI']
              # lower detection limit (expression below this 
              # will be treated as zero.
              lowerDetectionLimit: [0]
              # expressionFamily: for synthetic data use uninormal
              #  for mRNA count data use negbinomial.size()
              expressionFamily: ['uninormal']
              ### optional but recommended parameters
              # log transform expression values or not
              log: [False]
              # ignore pseudotime values (and use experimental
              # time points instead), recommended True for synthetic data
              # False for real mRNA data
              ignorePT: [True]