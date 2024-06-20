#This is the base image upon which necessary libraries are installed
FROM python:3.8

#Additional information
LABEL Maintainer="Yiqi Su <yiqisu@vt.edu>"

#Set main user as root to avoid permission issues
USER root

#Sets current working directory
WORKDIR /

#Copy the main Python script which will perform necessary GRN computations
COPY run_scSGL.py /

#Copy the original scSGL repo files stored in SCSGL folder
COPY scSGL /scSGL

#Installl time command to compute time taken
RUN apt-get update && apt-get install -y r-base time

#Install the requirments and install R to conda environment
RUN pip install -r /scSGL/requirements.txt

#Install pcaPP to use zero inflated Kendall tau as a kernel
RUN Rscript -e "install.packages('pcaPP')"

#Make a directory to mount the folder containing input files
RUN mkdir data/
