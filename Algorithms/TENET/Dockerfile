FROM python:3

LABEL Maintainer="Tim Wilson <timwilson@vt.edu>"

USER root

WORKDIR /

RUN apt-get update && apt-get install -y time && \ 
    apt-get install wget -y &&\
    apt-get install openmpi-bin -y


    # testing binary version  
    # wget https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.3.tar.bz2 && \
    # tar xf openmpi-5.0.3.tar.bz2 && \
    # cd openmpi-5.0.3 && \
    # ./configure --prefix=/configure
    

RUN pip install jpype1 statsmodels numpy

COPY TENET /TENET

RUN mkdir /data

WORKDIR /TENET