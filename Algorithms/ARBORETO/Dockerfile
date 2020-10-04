
FROM continuumio/anaconda3:2018.12

LABEL Maintainer="Aditya Pratapa <adyprat@vt.edu>"

USER root

RUN apt-get update

RUN conda install -y -c bioconda/label/cf201901 arboreto=0.1.5 pandas=0.24.0

COPY runArboreto.py /

RUN mkdir data/

RUN apt-get install time
