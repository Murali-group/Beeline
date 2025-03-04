FROM r-base:4.3.1

RUN apt-get update && apt-get install -y

USER root

WORKDIR /

# RUN R -e "install.packages(c( 'dplyr', 'optparse', 'purrr', 'rsample', 'igraph', 'parallel','tidyr'), repos='https://cloud.r-project.org/')"
# RUN R -e "install.packages(c('tidyverse', 'dplyr', 'optparse', 'purrr', 'rsample', 'igraph', 'parallel', 'splines', 'L0Learn','tidyr', 'BiocManager'), repos='https://cloud.r-project.org/')"
# RUN R -e "BiocManager::install('GENIE3')"
# EXPOSE 8080

RUN mkdir data/

COPY runMICA.R /
COPY spearman.R /
COPY l0l2.R /
COPY genie3.R /
COPY mutual_info.R /



LABEL maintainer="Vishnu Madhav <vm.ibab@gmail.com>"
LABEL version="1.0"
LABEL description="Mutual Information Chromatin Accessibility : Gene Regulatory network building"
