FROM r-base:3.5.3

LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

RUN apt-get update && apt-get install -y git

RUN apt-get install time

RUN R -e "install.packages(c('Rcpp', 'foreach', 'doParallel', 'poweRlaw', 'diptest', 'igraph', 'infotheo', 'entropy'))"

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/BTR_1.2.4.tar.gz', type = 'source', dependencies = TRUE)"

RUN mkdir data/

WORKDIR data/

COPY runBTR.R /data/runBTR.R



