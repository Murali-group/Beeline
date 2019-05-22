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

RUN mkdir runGRISLI

COPY runGRISLI/ /runGRISLI/

WORKDIR runGRISLI/

ENV LD_LIBRARY_PATH /opt/mcr/v96/runtime/glnxa64:/opt/mcr/v96/bin/glnxa64

RUN mkdir -p /root/.mcrCache9.6/GRISLI0/GRISLI/

RUN cp -r /runGRISLI/spams-matlab-v2.6/ /root/.mcrCache9.6/GRISLI0/GRISLI/

RUN mkdir data/

RUN apt-get update

RUN apt-get install -y libgomp1 --fix-missing

RUN apt-get install time