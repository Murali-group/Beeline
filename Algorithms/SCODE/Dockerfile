FROM r-base

LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

RUN apt-get update && apt-get install -y git

RUN apt-get install -y ruby

RUN git clone https://github.com/hmatsu1226/SCODE

WORKDIR SCODE/

RUN R -e "install.packages('MASS')"

RUN ruby run_R.rb data/exp_train.txt data/time_train.txt out 100 4 356 100 2


