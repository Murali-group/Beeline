# FROM r-base:4.2.2
FROM r-base:4.5.1

LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/LEAP/LEAP_0.2.tar.gz', type = 'source')"

COPY runLeap.R /

RUN mkdir data/

RUN apt-get update && apt-get install -y time
