FROM r-base

LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

RUN R -e "install.packages('glmnet')"
RUN R -e "install.packages('kSamples')"
RUN R -e "install.packages('ppcor')"
RUN R -e "install.packages('pracma')"
RUN R -e "install.packages('R.matlab')"
RUN R -e "install.packages('cvTools')"

RUN ls

COPY SINCERITIES.zip /


RUN unzip SINCERITIES.zip -d SINCERITIES

WORKDIR SINCERITIES/

RUN mkdir data/

RUN apt-get update && apt-get install time