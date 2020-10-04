FROM r-base

LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/glmnet/glmnet_3.0.tar.gz',type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/kSamples_1.2-9.tar.gz', type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/ppcor_1.1.tar.gz', type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/pracma_2.2.9.tar.gz',type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/R.matlab_3.6.2.tar.gz', type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/cvTools_0.3.2.tar.gz', type='source')"

RUN ls

COPY SINCERITIES.zip /


RUN unzip SINCERITIES.zip -d SINCERITIES

WORKDIR SINCERITIES/

RUN mkdir data/

RUN apt-get update && apt-get install time
