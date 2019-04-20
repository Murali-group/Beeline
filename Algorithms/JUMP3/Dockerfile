FROM ubuntu:18.04

LABEL Maintainer = "Aditya Pratapa <adyprat@vt.edu>"


RUN apt-get -qq update && apt-get -qq install -y unzip xorg wget curl libstdc++6

RUN mkdir /mcr-install && \
    mkdir /opt/mcr && \
    cd /mcr-install && \
    wget -q http://ssd.mathworks.com/supportfiles/downloads/R2019a/Release/0/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2019a_glnxa64.zip && \
    cd /mcr-install && \
    unzip MATLAB_Runtime_R2019a_glnxa64.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf mcr-install

RUN mkdir JUMP3/

COPY J3p/ /JUMP3/

WORKDIR JUMP3/

ENV LD_LIBRARY_PATH /opt/mcr/v96/runtime/glnxa64:/opt/mcr/v96/bin/glnxa64

RUN mkdir data/

RUN apt-get install time