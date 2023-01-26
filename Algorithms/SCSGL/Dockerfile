FROM python:3.8  #This is the base image upon which necessary libraries are installed

LABEL Maintainer="Yiqi Su <yiqisu@vt.edu>" #Additional information

USER root #Set main user as root to avoid permission issues

WORKDIR / #Sets current working directory

COPY run_scSGL.py /  #Copy the main Python script which will perform necessary GRN computations

COPY scSGL /scSGL  #Copy the original scSGL repo files stored in SCSGL folder 

RUN apt-get update && apt-get install -y r-base time  #Installl time command to compute time taken

RUN pip install -r /scSGL/requirements.txt #Install the requirments and install R to conda environment

RUN Rscript -e "install.packages('pcaPP')" #Install pcaPP to use zero inflated Kendall tau as a kernel

RUN mkdir data/ #Make a directory to mount the folder containing input files
