#!/bin/bash

# Building docker for the different algorithms 
echo "This may take a while..."

BASEDIR=$(pwd)

# You may remove the -q flag if you want to see the docker build status
cd $BASEDIR/Algorithms/ARBORETO
docker build -q -t arboreto:base .
if [[ "$(docker images -q arboreto:base 2> /dev/null)" != "" ]]; then
    echo "Docker container for ARBORETO is built and tagged as arboreto:base"
else
    echo "Oops! Unable to build Docker container for ARBORETO"
fi


cd $BASEDIR/Algorithms/GRISLI/
docker build -q -t grisli:base .
if [[ "$(docker images -q grisli:base 2> /dev/null)" != "" ]]; then
    echo "Docker container for GRISLI is built and tagged as grisli:base"
else
    echo "Oops! Unable to build Docker container for GRISLI"
fi


cd $BASEDIR/Algorithms/GRNVBEM/
docker build -q -t grnvbem:base .
if ([[ "$(docker images -q grnvbem:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for GRNVBEM is built and tagged as  grnvbem:base"
else
    echo "Oops! Unable to build Docker container for GRNVBEM"
fi


cd $BASEDIR/Algorithms/JUMP3/
docker build -q -t jump3:base .
if ([[ "$(docker images -q jump3:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for JUMP3 is built and tagged as  jump3:base"
else
    echo "Oops! Unable to build Docker container for JUMP3"
fi


cd $BASEDIR/Algorithms/LEAP/
docker build -q -t leap:base .
if ([[ "$(docker images -q leap:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for LEAP is built and tagged as  leap:base"
else
    echo "Oops! Unable to build Docker container for LEAP"
fi


cd $BASEDIR/Algorithms/PIDC/
docker build -q -t pidc:base .
if ([[ "$(docker images -q pidc:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for PIDC is built and tagged as pidc:base"
else
    echo "Oops! Unable to build Docker container for PIDC"
fi


cd $BASEDIR/Algorithms/PNI/
docker build -q -t pni:base .
if ([[ "$(docker images -q pni:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for PNI is built and tagged as pni:base"
else
    echo "Oops! Unable to build Docker container for PNI"
fi


cd $BASEDIR/Algorithms/PPCOR/
docker build -q -t ppcor:base .
if ([[ "$(docker images -q ppcor:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for PPCOR is built and tagged as ppcor:base"
else
    echo "Oops! Unable to build Docker container for PPCOR"
fi


cd $BASEDIR/Algorithms/SINGE/
docker build -q -t singe:base .
if ([[ "$(docker images -q singe:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for SINGE is built and tagged as singe:base"
else
    echo "Oops! Unable to build Docker container for SINGE"
fi


cd $BASEDIR/Algorithms/SCNS/
docker build -q -t scns:base .
if ([[ "$(docker images -q scns:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for SCNS is built and tagged as scns:base"
else
    echo "Oops! Unable to build Docker container for SCNS"
fi


cd $BASEDIR/Algorithms/SCODE/
docker build -q -t scode:base .
if ([[ "$(docker images -q scode:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for SCODE is built and tagged as scode:base"
else
    echo "Oops! Unable to build Docker container for SCODE"
fi


cd $BASEDIR/Algorithms/SCRIBE/
docker build -q -t scribe:base .
if [[ "$(docker images -q scribe:base 2> /dev/null)" != "" ]]; then
    echo "Docker container for SCRIBE is built and tagged as scribe:base"
else
    echo "Oops! Unable to build Docker container for SCRIBE"
fi


cd $BASEDIR/Algorithms/SINCERITIES/
docker build -q -t sincerities:base .
if ([[ "$(docker images -q sincerities:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for SINCERITIES is built and tagged as sincerities:base"
else
    echo "Oops! Unable to build Docker container for SINCERITIES"
fi


cd $BASEDIR
