FROM r-base

LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/versions_0.3.tar.gz', type='source')"
RUN R -e "require(versions); install.versions('glmnet', version='2.0-13')"
RUN R -e "require(versions); install.versions('kSamples', version='1.2-9')"	
RUN R -e "require(versions); install.versions('ppcor', version='1.1')"	
RUN R -e "require(versions); install.versions('pracma', version='2.2.9')"	
RUN R -e "require(versions); install.versions('R.matlab', version='3.6.2')"	
RUN R -e "require(versions); install.versions('cvTools', version='0.3.2')"	


RUN ls

COPY SINCERITIES.zip /


RUN unzip SINCERITIES.zip -d SINCERITIES

WORKDIR SINCERITIES/

RUN mkdir data/

RUN apt-get update && apt-get install time
