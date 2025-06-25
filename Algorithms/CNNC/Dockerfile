FROM nvidia/cuda:12.6.3-cudnn-devel-ubuntu20.04

LABEL Maintainer="Tim Wilson <timwilson@vt.edu>"
USER root
WORKDIR /
SHELL ["/bin/bash", "-c"]

RUN apt-get update
RUN apt-get install -y time wget

RUN mkdir -p /miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda3/miniconda.sh && \
    bash /miniconda3/miniconda.sh -b -u -p /miniconda3 && \
    rm /miniconda3/miniconda.sh
ENV PATH /miniconda3/bin:$PATH
ENV CONDA_OVERRIDE_CUDA="11.8"
ADD environment.yml /tmp/environment.yml
RUN conda config --set channel_priority flexible
RUN conda env create -f /tmp/environment.yml --debug
RUN echo "source activate CNNC" > ~/.bashrc
ENV PATH /miniconda3/envs/CNNC/bin:$PATH
RUN mkdir /data
ENV KERAS_BACKEND=torch
ENV NVIDIA_DRIVER_CAPABILITIES=all
RUN conda list -n CNNC > /list.txt
# ENV CUDA_VISIBLE_DEVICES = "0"
COPY cnncTrainModelKFold.py /
