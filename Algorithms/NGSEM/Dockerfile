FROM r-base:4.4.1

LABEL maintainer = "Tim Wilson <timwilson@vt.edu>"

USER root

WORKDIR /

# RUN apt-get update && apt-get install --allow-unauthenticated -y gnupg2
# RUN apt-key net-update && apt-get update && apt-get install --allow-unauthenticated -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev time

RUN apt-get update --allow-unauthenticated && apt-get install -y --allow-unauthenticated usrmerge libcurl4-gnutls-dev libxml2-dev libssl-dev build-essential
#RUN R -e 'install.packages("pkgdown")' && R -e 'install.packages("devtools", dependencies=TRUE)' && R -e 'install_github("jiayingzhao/NGSEM")'

COPY NGSEM /
RUN R -e "install.packages('NGSEM', repos = NULL, type='source')"
COPY runNGSEM.R /
RUN mkdir data/