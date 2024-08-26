FROM anibali/pytorch:1.7.0-cuda11.0

USER root

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Moscow

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
        tzdata \
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		fonts-texgyre \
        software-properties-common\
		sudo \
		make \
		build-essential \
	&& apt-get clean && rm -rf /var/lib/apt/lists/*

ENV R_BASE_VERSION 4.1.2

RUN sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
	&& apt-get update \
        && apt-get install -y --no-install-recommends \
                libopenblas0-pthread \
		littler \
                r-cran-littler \
		r-base=${R_BASE_VERSION}-* \
		r-base-dev=${R_BASE_VERSION}-* \
                r-base-core=${R_BASE_VERSION}-* \
		r-recommended=${R_BASE_VERSION}-* \
	&& ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/installBioc.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installDeps.r /usr/local/bin/installDeps.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
	&& install.r docopt \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*

# rocker/r-base end
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE

# Install Seurat's system dependencies
RUN apt-get update \
	&& apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libboost-all-dev \
    libxml2-dev \
    openjdk-8-jdk \
    python3-dev \
    python3-pip \
    git \
    libfftw3-dev \
    libgsl-dev \
	llvm-10


# Install UMAP
RUN LLVM_CONFIG=/usr/lib/llvm-10/bin/llvm-config pip3 install llvmlite \
	&& pip3 install numpy==1.21 \
	&& pip3 install umap-learn

# Install FIt-SNE
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git
RUN g++ -std=c++11 -O3 FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp  -o FIt-SNE/bin/fast_tsne -pthread -lfftw3 -lm

# Install bioconductor dependencies & suggests
RUN R --no-echo --no-restore --no-save -e "install.packages('BiocManager')" \
	&& R --no-echo --no-restore --no-save -e "BiocManager::install(c('multtest', 'S4Vectors', 'SummarizedExperiment', 'SingleCellExperiment', 'MAST', 'DESeq2', 'BiocGenerics', 'GenomicRanges', 'IRanges', 'rtracklayer', 'monocle', 'Biobase', 'limma'))"

# Install CRAN suggests
RUN R --no-echo --no-restore --no-save -e "install.packages(c('VGAM', 'R.utils', 'metap', 'Rfast2', 'ape', 'enrichR', 'mixtools'))"

# Install hdf5r
RUN R --no-echo --no-restore --no-save -e "install.packages('hdf5r')"

# Install Seurat
RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')" \
	&& R --no-echo --no-restore --no-save -e "install.packages('Seurat')"

# Install SeuratDisk
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('mojaveazure/seurat-disk')"

## seurat end

WORKDIR /tmp

RUN NCPUS=${NCPUS:-1}
RUN set -e
RUN apt-get -q  update && apt-get -y --no-install-recommends install \
	libxml2-dev \
	libcairo2-dev \
	libgit2-dev \
	default-libmysqlclient-dev \
	libpq-dev \
	libsasl2-dev \
	libsqlite3-dev \
	libssh2-1-dev \
	unixodbc-dev \
	python \
	python3-pip \
	libbz2-dev \
	liblzma-dev \
	libsodium-dev \
	libhiredis-dev

###############
# HTSlib 1.11.0#
# Folk from https://github.com/chrisamiller/docker-r-seurat/blob/master/Dockerfile
###############
#ENV HTSLIB_INSTALL_DIR=/opt/htslib

RUN wget --no-check-certificate https://github.com/samtools/htslib/archive/1.11.0.zip && \
	unzip 1.11.0.zip && \
	rm 1.11.0.zip && \
	cd /tmp/htslib-1.11.0 && \
	#./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
	make && \
	make install && \ 
	cp -R * /usr/lib/

# Install Bioconductor dependencies and Install Bioconductor databases and # Install CRAN dependencies
RUN R -e 'BiocManager::install(c("GenomicAlignments", "ggbio", "biovizBase", "fgsea", "ComplexHeatmap", "karyoploteR", "MAGeCKFlute"),force = TRUE)' \
	&& R -e 'BiocManager::install(c("JASPAR2020", "GO.db", "EnsDb.Hsapiens.v86","EnsDb.Mmusculus.v79", "org.Hs.eg.db", "org.Mm.eg.db","scater","bluster","dsb","tidyverse"),force = TRUE)' \
	&& install2.r --error --skipinstalled \
	Polychrome \
	qs \
	plumber \
	vroom \
	lintr \
	gert \
	Signac \ 
	logger \
	tictoc \
	msigdbr \
	Gmisc \ 
	rematch \
	readr \
	openxlsx \
	readxl \
	sp \
	statmod \
	statmod\
	nloptr\
	minqa\
	lme4\
	rio\
	maptools\
	pbkrtest\
	carData\
	car\
	corrplot\
	broom\
	rstatix\
	polynom\
	ggsignif\
	ggsci\
	ggpubr \
	spatstat \
	rlist \
	redux \
	devtools \
	enrichR \
	NMF \
	ggalluvial  \
	svglite \
	expm \
	sna \
	gg.gap \
	&& installGithub.r -d FALSE -u FALSE\
	immunogenomics/presto \ 
	liulab-dfci/MAESTRO \
	sqjin/CellChat 

# Socket-io 
RUN pip3 install python-socketio[client]==4.6.1 \
	&& rm -rf /tmp/*  \
	&& rm -rf /var/lib/apt/lists/*

# CMD ["Rscript", "-e", "installed.packages()"]

# r base end

# Python
RUN conda update -y conda \
	&& conda install -y -c anaconda pip \
	#&& pip install -U pip \
	&& conda config --add channels conda-forge \
    && conda install numpy numpy_groupies scipy matplotlib pandas seaborn scikit-learn notebook dash plotly black bokeh h5py click jupyter jupyterlab pytables \
	&& pip install -U --no-cache-dir jupyterthemes jupyter_contrib_nbextensions python-igraph umap-learn numba Cython transformers pyreadr dill redis\
	&& jupyter notebook --generate-config \
	&& jupyter lab clean

# Install Pytorch Geometric, Velocity and kneed
RUN pip install --no-cache-dir torch-scatter==2.0.7 -f https://data.pyg.org/whl/torch-1.7.0+cu110.html \
    && pip install torch-sparse==0.6.9 -f https://data.pyg.org/whl/torch-1.7.0+cu110.html \
    && pip install torch-cluster -f https://data.pyg.org/whl/torch-1.7.0+cu110.html \
    && pip install torch-spline-conv -f https://data.pyg.org/whl/torch-1.7.0+cu110.html \
    && pip install torch-geometric \
    && pip install -U --no-cache-dir \
	velocyto \
	scvelo \
    kneed


# socket-io must be v4.6
RUN pip install lisa2==2.2.5 \
	&& pip install -U --no-cache-dir python-socketio==4.6.1 \
	&& conda clean --all -f -y 

# Add Tini
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini

EXPOSE 8888
# ENTRYPOINT ["/tini", "--"]
# CMD ["/bin/bash"]