FROM r-base:4.2.2

LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/versions_0.3.tar.gz', type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/glmnet/glmnet_2.0-13.tar.gz', type='source')"
RUN R -e "install.packages('https://packagemanager.posit.co/cran/latest/src/contrib/Archive/kSamples/kSamples_1.2-9.tar.gz', type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/ppcor_1.1.tar.gz', type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/pracma/pracma_2.2.9.tar.gz', type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/R.matlab/R.matlab_3.6.2.tar.gz', type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/cvTools/cvTools_0.3.2.tar.gz', type='source')"


RUN ls

COPY SINCERITIES.zip /


RUN unzip SINCERITIES.zip -d SINCERITIES

WORKDIR SINCERITIES/

RUN mkdir data/

RUN apt-get update && apt-get install -y time
