FROM r-base:4.2.0

LABEL maintainer = "Daniel Osorio <dcosorioh@tamu.edu>"

USER root

WORKDIR /

# RUN R -e "install.packages('https://cran.r-project.org/src/contrib/remotes_2.4.2.tar.gz', type = 'source')"

RUN R -e "install.packages('remotes')" 

RUN R -e "library(remotes)"

RUN R -e "remotes::install_cran(pkgs = 'scTenifoldNet', quiet = TRUE)"

RUN R -e "remotes::install_cran(pkgs = 'reshape2', quiet = TRUE)"

# RUN R -e "install.packages('scTenifoldNet')" 

# RUN R -e "install.packages('reshape2')"

RUN R -e "library(scTenifoldNet)"

RUN R -e "library(reshape2)"

COPY runSCTENIFOLDNET.R /

RUN mkdir data/

RUN apt-get update && apt-get install -y time

