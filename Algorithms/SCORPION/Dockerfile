FROM r-base:4.2.0

LABEL maintainer = "Daniel Osorio <daniel.osoriohurtado@childrens.harvard.edu>"

USER root

WORKDIR /

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/remotes_2.4.2.tar.gz', type = 'source')"

RUN R -e "install.packages('reshape2')"

# RUN R -e "remotes::install_github('kuijjerlab/SCORPION')"

RUN R -e "install.packages('SCORPION')"

RUN R -e "library(reshape2)"

RUN R -e "library(SCORPION)"

COPY runSCORPION.R /

RUN mkdir data/

RUN apt-get update && apt-get install -y time
