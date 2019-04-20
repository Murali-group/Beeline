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

RUN apt-get install time