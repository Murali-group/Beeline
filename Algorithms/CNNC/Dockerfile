FROM continuumio/miniconda3

LABEL Maintainer="Tim Wilson <timwilson@vt.edu>"
USER root
WORKDIR /
# RUN apt-get update
# RUN apt-get install -y curl gpg
# RUN curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg \
#   && curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
#     sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
#     tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
RUN apt-get update
RUN apt-get install -y time
# RUN apt-get install -y nvidia-container-toolkit

ADD environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml
RUN echo "source activate CNNC" > ~/.bashrc
ENV PATH /opt/conda/envs/CNNC/bin:$PATH
RUN mkdir /data
ENV KERAS_BACKEND=theano
RUN mkdir /files/
COPY cnncTrainModelKFold.py /
# COPY cnncDataReprod.py /
# COPY NEPDF_data /NEPDF_data
# can remove these after reprod
# RUN conda install pandas
# RUN conda install numpy
# RUN conda install matplotlib 
# RUN conda install scikit-learn 
# RUN conda install scipy 
# RUN conda install theano
# RUN conda install keras
