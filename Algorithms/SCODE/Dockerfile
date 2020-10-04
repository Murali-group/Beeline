FROM r-base:3.5.3

LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /


RUN apt-get update && apt-get install -y  gcc-8-base 

RUN apt-get install -y git

RUN apt-get install -y ruby

RUN git clone https://github.com/hmatsu1226/SCODE

WORKDIR SCODE/

RUN git checkout a0512f8ec29aac188c9c27a8e89ddd2464e6d84d

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/MASS_7.3-51.3.tar.gz', type = 'source')"

# RUN ruby run_R.rb data/exp_train.txt data/time_train.txt out 100 4 356 100 2



RUN apt-get install time
