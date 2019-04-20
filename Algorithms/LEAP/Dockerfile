FROM r-base:3.5.3

LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/LEAP_0.2.tar.gz', type = 'source')"

COPY runLeap.R /

RUN mkdir data/

RUN apt-get update && apt-get install time