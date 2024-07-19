FROM python:3

LABEL Maintainer="Tim Wilson <timwilson@vt.edu>"

USER root

WORKDIR /

RUN apt-get update && apt-get install -y time && \ 
    apt-get install wget -y &&\
    # install openmpi
    wget https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.3.tar.bz2 && \
    tar xf openmpi-5.0.3.tar.bz2 && \
    cd openmpi-5.0.3 && \
    mkdir build && \
    cd build && \
    ../configure && \
    make -j 4 && \
    apt-get install default-jre -y && \
    apt-get install default-jdk -y
    

RUN pip install jpype1 statsmodels numpy    

COPY TENET /TENET

RUN mkdir /data

WORKDIR /TENET