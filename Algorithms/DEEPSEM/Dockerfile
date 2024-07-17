FROM python:3.7

USER root

WORKDIR /

COPY /DEEPSEM /

RUN apt-get update && apt-get install time

RUN pip install wheel


RUN pip install torch==1.2.0 torchvision==0.4.0

RUN pip install 'scanpy[leiden]'

RUN pip install -r requirements.txt
    
RUN mkdir data/