#!/bin/bash

### Create the tar gzipped files with the C++ sources

# file name with date and version
DA=$(``date +%F)
echo $DA
VERSION=`cat swig/Version`
WDIR="spams-cpp-v$VERSION"
WFILE="spams-cpp-v$VERSION-$da.tar.gz"

# create the tar gzipped file
git archive --format=tar.gz --prefix=$WDIR --output=$WFILE HEAD
