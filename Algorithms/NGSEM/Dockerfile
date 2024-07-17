FROM r-base:4.2.2

LABEL maintainer = "Tim Wilson <timwilson@vt.edu>"

USER root

WORKDIR /

RUN apt-get update --allow-unauthenticated && apt-get install -y --allow-unauthenticated time r-cran-devtools

RUN R -e "install.packages('MASS')"
RUN R -e "install.packages('foreach')"
RUN R -e "install.packages('knitr')"
RUN R -e "install.packages('rmarkdown')"

RUN mkdir ng-sem/

WORKDIR ng-sem/

COPY NGSEM ng-sem/

COPY runNGSEM.R /ng-sem/

RUN mkdir data/