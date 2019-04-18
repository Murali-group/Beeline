FROM bioconductor/release_base2:R3.5.3_Bioc3.8

LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

Run R -e "BiocManager::install('devtools',version=3.8)"

RUN apt-get update && apt-get install -y r-cran-rgl

RUN R -e "BiocManager::install('HiveR',ref = 3.8)"


