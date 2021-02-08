# MATLAB R2018a in a Debian environment
FROM amarburg/matlab-runtime:R2018a

USER root

RUN apt-get update && \
    apt-get -y install libxt6 bzip2 time octave && \
    rm -rf /var/lib/apt/lists/*

# Install Miniconda3 following https://hub.docker.com/r/continuumio/miniconda3/dockerfile
# Python is only needed for testing SINGE and could be removed from the base
# to reduce the image size
ENV PATH /opt/conda/bin:$PATH
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean --all && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    conda init bash

ENV SINGE_ROOT /usr/local/SINGE
WORKDIR /usr/local/SINGE
ENTRYPOINT ["/usr/local/SINGE/SINGE.sh", "/usr/local/MATLAB/MATLAB_Runtime/v94"]

# Install conda test environment
COPY tests/environment.yml tests/
RUN conda env create -f tests/environment.yml

# Copy wrapper scripts for compiled MATLAB executables
COPY SINGE.sh .
COPY run_SINGE_GLG_Test.sh .
COPY run_SINGE_Aggregate.sh .
COPY tests/run_SINGE_Test.sh tests/

ENV SINGE_RUNNING_IN_DOCKER 1

# Download the compiled SINGE executables from the stable release
# md5sum of v0.4.1 SINGE_GLG_Test is c50ec7bc13e287eca340c9d19d8bc27d/
RUN tag=v0.4.1 && \
    wget --quiet https://github.com/gitter-lab/SINGE/releases/download/$tag/SINGE_Test && \
    wget --quiet https://github.com/gitter-lab/SINGE/releases/download/$tag/SINGE_GLG_Test && \
    wget --quiet https://github.com/gitter-lab/SINGE/releases/download/$tag/SINGE_Aggregate && \
    wget --quiet https://github.com/gitter-lab/SINGE/releases/download/$tag/code.md5 && \
    chmod u+x SINGE_* && \
    mv SINGE_Test tests/SINGE_Test

# Download an intermediate version of the compiled SINGE executables for testing
# Download the md5sums of the source .m files and binaries
#RUN md5=c50ec7bc13e287eca340c9d19d8bc27d && \
#    wget --quiet https://www.biostat.wisc.edu/~gitter/tmp/$md5/SINGE_Test && \
#    wget --quiet https://www.biostat.wisc.edu/~gitter/tmp/$md5/SINGE_GLG_Test && \
#    wget --quiet https://www.biostat.wisc.edu/~gitter/tmp/$md5/SINGE_Aggregate && \
#    wget --quiet https://www.biostat.wisc.edu/~gitter/tmp/$md5/code.md5 && \
#    chmod u+x SINGE_* && \
#    mv SINGE_Test tests/SINGE_Test
