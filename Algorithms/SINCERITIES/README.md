*This README.md file was generated on 1/20/2023 by Yiqi Su (yiqisu@vt.edu)*
<!-- remove all comments (like this) before final save  -->

# SINCERITIES: inferring gene regulatory networks from time-stamped single cell transcriptional expression profiles

This is the instruction on how to integrate SINCERITIES ([[Paper](https://doi.org/10.1093/bioinformatics/btx575)] [[GitHub](https://github.com/CABSEL/SINCERITIES)]) to BEELINE. Please follow the following steps:

1. **Create SINCERITIES folder:** Create a folder called SINCERITIES under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.

2. **Download SINCERITIES.zip file:** In the SINCERITIES folder, download a zip file SINCERITIES-R v2.0.zip from GitHub (https://github.com/CABSEL/SINCERITIES/blob/master/SINCERITIES-R_v2.0.zip) and keep the SINCERITIES function folder, MAIN.R script and ExpressionData.csv file to learn graphs from target datasets. 
   

3. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. 

        FROM r-base
        
        LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"
        
        USER root
        
        WORKDIR /
        
        RUN R -e "install.packages('https://cran.r-project.org/src/contrib/versions_0.3.tar.gz', type='source')"
        RUN R -e "require(versions); install.versions('glmnet', version='2.0-13')"
        RUN R -e "require(versions); install.versions('kSamples', version='1.2-9')"	
        RUN R -e "require(versions); install.versions('ppcor', version='1.1')"	
        RUN R -e "require(versions); install.versions('pracma', version='2.2.9')"	
        RUN R -e "require(versions); install.versions('R.matlab', version='3.6.2')"	
        RUN R -e "require(versions); install.versions('cvTools', version='0.3.2')"	
        
        RUN ls
        
        COPY SINCERITIES.zip /
        
        
        RUN unzip SINCERITIES.zip -d SINCERITIES
        
        WORKDIR SINCERITIES/
        
        RUN mkdir data/
        
        RUN apt-get update && apt-get install -y time


The Dockerfile will run the SINCERITIES algorithm within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for SINCERITIES.

        cd $BASEDIR/Algorithms/SINCERITIES/
        docker build -q -t sincerities:base .
        if ([[ "$(docker images -q sincerities:base 2> /dev/null)" != "" ]]); then
            echo "Do
            cker container for SINCERITIES is built and tagged as sincerities:base"
        else
            echo "Oops! Unable to build Docker container for SINCERITIES"
        fi

5. **Create sinceritiesRunner.py script:** After buliding the Docker image, create a Python script called sinceritiesRunner.py in Beeline/BLRun folder to setup BLRun objects to read inputs and run SINCERITIES inside the Docker image, and also parse the output for evaluation. Specifically, the sinceritiesRunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the input data files (i.e., expression data and pseudotime data), and processes them into the format required by the SINCERITIES. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data files (i.e., expression data and pseudotime data). It also specifies where the outputs are written. The docker container run SINCERITIES when the parameters are passed. 
   - ``parseOutput()`` : This function reads the sincerities-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

6. **Add SINCERITIES to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to SINCERITIES. 

    - add "import BLRun.sinceritiesRunner as SINCERITIES"
    - add "'SINCERITIES':SINCERITIES.generateInputs" to InputMapper
    - add "'SINCERITIES':SINCERITIES.run" to AlgorithmMapper
    - add "'SINCERITIES':SINCERITIES.parseOutput" to OutputParser

7. **Add SINCERITIES to config.yaml:** The final step is to add the new algorithm SINCERITIES and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.

        - name: "SINCERITIES"
          params: 
              should_run: [True]
              nBins: [10]