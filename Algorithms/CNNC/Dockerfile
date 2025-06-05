FROM continuumio/miniconda3

LABEL Maintainer="Tim Wilson <timwilson@vt.edu>"
USER root
WORKDIR /

ADD environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml
RUN echo "source activate CNNC" > ~/.bashrc
ENV PATH /opt/conda/envs/CNNC/bin:$PATH
RUN mkdir /data
RUN apt-get update
RUN apt-get install time
ENV KERAS_BACKEND=theano
RUN mkdir /files/
COPY cnncTrainModelKFold.py /
COPY cnncDataReprod.py /
COPY NEPDF_data /NEPDF_data
# can remove these after reprod
# RUN conda install pandas
# RUN conda install numpy
# RUN conda install matplotlib 
# RUN conda install scikit-learn 
# RUN conda install scipy 
# RUN conda install theano
# RUN conda install keras
