#!/bin/bash

# Building docker for the different algorithms 
echo "This may take a while..."

BASEDIR=$(pwd)
DOCKER_IMAGE_DIR=${BASEDIR}/images_docker
mkdir -p ${DOCKER_IMAGE_DIR}
SIF_DIR=${BASEDIR}/images_singularity
mkdir -p ${SIF_DIR}

#for ALG in ARBORETO GRISLI GRNVBEM JUMP3 LEAP PIDC PNI PPCOR SINGE SCNS SCODE SCRIBE SINCERITIES ; do
# for ALG in PPCOR SINGE SCNS SCODE SCRIBE SINCERITIES ; do
for ALG in SINGE ; do
    cd $BASEDIR/Algorithms/${ALG}
    docker build -q -t ${ALG,,}:base .
    echo "Docker container for ${ALG} is built and tagged as ${ALG,,}:base"
    DOCKER_ARCHIVE="${DOCKER_IMAGE_DIR}/${ALG}.tar"
    docker save ${ALG,,}:base -o ${DOCKER_ARCHIVE}
    echo "Docker container for ${ALG} is exported to ${DOCKER_ARCHIVE}"
    SIF_IMAGE="${SIF_DIR}/${ALG}.sif"
    singularity build ${SIF_IMAGE} docker-archive://${DOCKER_ARCHIVE}
    echo "Singularity image for ${ALG} is built and saved to ${SIF_IMAGE}"
done
cd $BASEDIR

cd $SIF_DIR
ln -s ARBORETO.sif GENIE3.sif
ln -s ARBORETO.sif GRNBOOST2.sif
cd $BASEDIR
