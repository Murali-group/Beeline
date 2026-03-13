FROM r-base:4.2.2

LABEL maintainer="Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

RUN R -e "install.packages('ppcor', repos='https://cloud.r-project.org', dependencies=TRUE)"

COPY runPPCOR.R /

RUN mkdir data/

RUN apt-get update -o Acquire::AllowInsecureRepositories=true \
    && apt-get install -y --allow-unauthenticated debian-archive-keyring \
    && apt-get update \
    && apt-get install -y time
