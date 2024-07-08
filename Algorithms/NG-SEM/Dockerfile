FROM r-base

LABEL maintainer = "Tim Wilson <timwilson@vt.edu>"

USER root

WORKDIR /

RUN R -e 'install.packages("devtools")' && R -e 'install_github("jiayingzhao/NG-SEM")'

COPY runNG-SEM.R /
RUN mkdir data/

RUN apt-get update && apt-get install -y time