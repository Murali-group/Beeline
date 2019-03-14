FROM ubuntu:18.04

LABEL maintainer="Aditya Pratapa <adyprat@vt.edu>"
USER root

RUN apt-get update && apt-get install -y apt-utils \
    	git\
    	wget \
    	tar \
    	make \
    	python \
	python-pip \
	python-dev \
	gcc \
