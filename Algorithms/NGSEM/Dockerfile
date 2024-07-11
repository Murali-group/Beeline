FROM r-base:4.2.2

LABEL maintainer = "Tim Wilson <timwilson@vt.edu>"

USER root

WORKDIR /

# RUN apt-get update && apt-get install --allow-unauthenticated -y gnupg2
# RUN apt-key net-update && apt-get update && apt-get install --allow-unauthenticated -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev time

RUN apt-get update --allow-unauthenticated && apt-get install -y --allow-unauthenticated time r-cran-devtools

RUN R -e "library(devtools); \ 
          install_github('jiayingzhao/NG-SEM')"
# COPY NGSEM /
#RUN R -e "install.packages('/NGSEM', repos = NULL, type='source')"
#RUN R -e "installed.packages()"
COPY runNGSEM.R /

RUN mkdir data/