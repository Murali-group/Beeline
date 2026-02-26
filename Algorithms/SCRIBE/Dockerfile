FROM bioconductor/bioconductor_docker:RELEASE_3_14

LABEL maintainer="Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

RUN R -e "install.packages('devtools')"

RUN apt-get update && apt-get install -y libhdf5-dev \
libxml2-dev \
libudunits2-dev \
imagemagick \
zlib1g-dev \
libfreetype6-dev

RUN DEBIAN_FRONTEND=noninteractive apt-get -y install xorg \
libx11-dev \
libglu1-mesa-dev

RUN R -e "install.packages('rgl')"


RUN R -e "BiocManager::install('HiveR', update = FALSE)"

RUN R -e "BiocManager::install(c('lattice','Matrix','irlba'), update = FALSE)"

RUN R -e "BiocManager::install(c('cluster','dplyr'), update = FALSE)"

# monocle 2.x — available in Bioconductor 3.14 as version 2.22.0

RUN R -e "BiocManager::install('monocle', update = FALSE)"

RUN git clone https://github.com/cole-trapnell-lab/RANNinf

WORKDIR /RANNinf/

RUN git checkout 9d4a3d781c8b74a01ec2018ca7137b22a0e583b6

WORKDIR /

RUN R CMD build RANNinf

RUN R -e "BiocManager::install('RcppArmadillo', update = FALSE)"

RUN R -e "install.packages('RANNinf_2.5.0.99.tar.gz', repo = NULL, type ='source')"

RUN git clone https://github.com/cole-trapnell-lab/Scribe

WORKDIR /Scribe/

RUN git checkout 4ba98500764adbce4a59be508d94b279bbfcfb31

RUN R -e "BiocManager::install(c('cowplot','lpSolveAPI'), update = FALSE)"

RUN R -e "install.packages('Scribe_0.1.tar.gz', repo = NULL, type ='source')"

WORKDIR /

RUN R -e "BiocManager::install(c('optparse'), update = FALSE)"

RUN mkdir data/

COPY runScribe.R /

RUN apt-get update && apt-get install -y time
