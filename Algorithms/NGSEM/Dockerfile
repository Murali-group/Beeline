FROM r-base:4.2.2

LABEL maintainer = "Tim Wilson <timwilson@vt.edu>"

USER root

WORKDIR /

# RUN apt-get update && apt-get install --allow-unauthenticated -y gnupg2
# RUN apt-key net-update && apt-get update && apt-get install --allow-unauthenticated -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev time

RUN apt-get update --allow-unauthenticated && apt-get install -y --allow-unauthenticated time r-cran-devtools git
COPY NG-SEM /
WORKDIR /NG-SEM
# RUN R CMD BUILD NG-SEM
# RUN R CMD INSTALL NG-SEM_0.1.tar.gz
RUN R -e "library(devtools); \
        build(); \
        install()"
# RUN R -e "library(devtools); \ 
#            install_github('jiayingzhao/NG-SEM')"


WORKDIR /
# RUN R -e "install.packages('https://github.com/jiayingzhao/NG-SEM/archive/master.tar.gz', type = 'source')"
#RUN R -e "install.packages('/NGSEM', repos = NULL, type='source')"
#RUN R -e "installed.packages()"
COPY runNGSEM.R /

RUN mkdir data/

# RUN mkdir R/

# RUN mv NG-SEM/R/* R/

# RUN mkdir man/

# RUN mv NG-SEM/man/* man/

# RUN mv NG-SEM/NAMESPACE /