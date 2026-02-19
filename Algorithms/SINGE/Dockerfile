# MATLAB R2018a in a Debian environment
FROM amarburg/matlab-runtime:R2018a

RUN echo "deb [trusted=yes] http://archive.debian.org/debian stretch main non-free contrib" > /etc/apt/sources.list && \
    echo "deb-src [trusted=yes] http://archive.debian.org/debian stretch main non-free contrib" >> /etc/apt/sources.list && \
    echo "deb [trusted=yes] http://archive.debian.org/debian-security stretch/updates main non-free contrib" >> /etc/apt/sources.list

# RUN gpg --keyserver keyserver.ubuntu.com --recv-keys AA8E81B4331F7F50 && \
#     gpg --export --armor AA8E81B4331F7F50 | sudo tee /etc/apt/trusted.gpg.d/AA8E81B4331F7F50.gpg

RUN cat /etc/apt/sources.list

RUN apt-get update && \
    apt-get -y install libxt6 bzip2 && \
    rm -rf /var/lib/apt/lists/*

# Install Miniconda3 following https://hub.docker.com/r/continuumio/miniconda3/dockerfile
# Python is only needed for testing SINGE and could be removed from the base
# to reduce the image size
ENV PATH=/opt/conda/bin:$PATH
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean --all && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    conda init bash

ENV SINGE_ROOT=/usr/local/SINGE
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

ENV SINGE_RUNNING_IN_DOCKER=1

# Download the compiled SINGE executables from the stable release
# md5sum of v0.5.0 SINGE_GLG_Test is 1834d44cb686cdfd9b7ee99398518115
RUN tag=v0.5.0 && \
   wget --quiet https://github.com/gitter-lab/SINGE/releases/download/$tag/SINGE_Test && \
   wget --quiet https://github.com/gitter-lab/SINGE/releases/download/$tag/SINGE_GLG_Test && \
   wget --quiet https://github.com/gitter-lab/SINGE/releases/download/$tag/SINGE_Aggregate && \
   wget --quiet https://github.com/gitter-lab/SINGE/releases/download/$tag/code.md5 && \
   chmod u+x SINGE_* && \
   mv SINGE_Test tests/SINGE_Test

# Download an intermediate version of the compiled SINGE executables for testing
# Download the md5sums of the source .m files and binaries
# RUN md5=e9bc31dfc6274f4874cc2bde00d85787 && \
#     wget --quiet --no-check-certificate https://www.biostat.wisc.edu/~gitter/tmp/$md5/SINGE_Test && \
#     wget --quiet --no-check-certificate https://www.biostat.wisc.edu/~gitter/tmp/$md5/SINGE_GLG_Test && \
#     wget --quiet --no-check-certificate https://www.biostat.wisc.edu/~gitter/tmp/$md5/SINGE_Aggregate && \
#     wget --quiet --no-check-certificate https://www.biostat.wisc.edu/~gitter/tmp/$md5/code.md5 && \
#     chmod u+x SINGE_* && \
#     mv SINGE_Test tests/SINGE_Test