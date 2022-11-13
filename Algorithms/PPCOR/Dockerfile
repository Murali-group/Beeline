FROM r-base:4.2.2

LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/ppcor_1.1.tar.gz', type = 'source')"

COPY runPPCOR.R /

RUN mkdir data/

RUN apt-get update && apt-get install -y time