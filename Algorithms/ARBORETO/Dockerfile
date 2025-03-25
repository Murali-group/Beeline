FROM continuumio/anaconda3:2024.02-1

LABEL Maintainer="Aditya Pratapa <adyprat@vt.edu>"

USER root

RUN apt-get update

RUN conda install -y -c bioconda arboreto pandas

COPY runArboreto.py /

RUN mkdir data/

RUN apt-get install time
