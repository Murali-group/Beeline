#!/bin/bash

set -e # abandon script on error
BASEDIR="$(dirname "$(readlink -f "$0")")" # set env variable for current directory

BUILD=false
HELP=false
VERBOSE_VALUE="-q "

show_help() {
  echo "Usage: $(basename "$0") [OPTIONS] [ARGUMENTS]"
  echo "This script creates docker containers for BEELINE."
  echo ""
  echo "Options:"
  echo "  -h, --help    Display this help message and exit."
  echo "  -b, --build   Instead of pulling images from docker hub, build them manually locally."
  echo "  -v, --verbose Enable verbose output."
  echo ""
  echo "Requirements:"
  echo "  docker (last version tested 28.5.1, build e180ab8)"
  echo "  conda (last version tested 25.9.0)"
  echo "  git"
  echo ""
  echo "Examples:"
  echo "  $(basename "$0")"
  echo "  $(basename "$0") -b --verbose"
}

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -b|--build)
      BUILD=true
      ;;
    -v|--verbose)
      VERBOSE_VALUE=""
      ;;
    -h|--help)
      HELP=true
      ;;
    *)
      echo "Unknown option: $1" >&2
      show_help
      exit 1
      ;;
  esac
  shift # Consume the current argument (flag or value)
done

if [[ "$HELP" = true ]]; then
    show_help
    exit 0
fi

if [[ "$BUILD" = true ]]; then

    # Building docker for the different algorithms 
    echo "This may take a while..."
    
    # You may remove the -q flag if you want to see the docker build status
    pushd $BASEDIR/Algorithms/ARBORETO
    docker build -t arboreto:base .
    if ([ $? = 0 ] && [ "$(docker images -q arboreto:base 2> /dev/null)" != "" ]); then
        echo "Docker container for ARBORETO is built and tagged as arboreto:base"
    elif [ "$(docker images -q arboreto:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at arboreto:base"
    else
        echo "Oops! Unable to build Docker container for ARBORETO"
    fi
    popd

    pushd $BASEDIR/Algorithms/GRISLI/
    docker build -t grisli:base .
    if ([ $? = 0 ] && [ "$(docker images -q grisli:base 2> /dev/null)" != "" ]); then
        echo "Docker container for GRISLI is built and tagged as grisli:base"
    elif [ "$(docker images -q grisli:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at grisli:base"
    else
        echo "Oops! Unable to build Docker container for GRISLI"
    fi
    popd

    pushd $BASEDIR/Algorithms/GRNVBEM/
    docker build -t grnvbem:base .
    if ([ $? = 0 ] && [[ "$(docker images -q grnvbem:base 2> /dev/null)" != "" ]]); then
        echo "Docker container for GRNVBEM is built and tagged as  grnvbem:base"
    elif [ "$(docker images -q grnvbem:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at grnvbem:base"
    else
        echo "Oops! Unable to build Docker container for GRNVBEM"
    fi
    popd

    pushd $BASEDIR/Algorithms/JUMP3/
    docker build -t jump3:base .
    if ([ $? = 0 ] && [[ "$(docker images -q jump3:base 2> /dev/null)" != "" ]]); then
        echo "Docker container for JUMP3 is built and tagged as  jump3:base"
    elif [ "$(docker images -q jump3:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at jump3:base"
    else
        echo "Oops! Unable to build Docker container for JUMP3"
    fi
    popd

    pushd $BASEDIR/Algorithms/LEAP/
    docker build --tag=leap:base .
    if ([ $? = 0 ] && [[ "$(docker images -q leap:base 2> /dev/null)" != "" ]]); then
        echo "Docker container for LEAP is built and tagged as  leap:base"
    elif [ "$(docker images -q leap:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at leap:base"
    else
        echo "Oops! Unable to build Docker container for LEAP"
    fi
    popd

    pushd $BASEDIR/Algorithms/PIDC/
    docker build -t pidc:base .
    if ([ $? = 0 ] && [[ "$(docker images -q pidc:base 2> /dev/null)" != "" ]]); then
        echo "Docker container for PIDC is built and tagged as pidc:base"
    elif [ "$(docker images -q pidc:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at pidc:base"
    else
        echo "Oops! Unable to build Docker container for PIDC"
    fi
    popd

    pushd $BASEDIR/Algorithms/PNI/
    docker build -t pni:base .
    if ([ $? = 0 ] && [[ "$(docker images -q pni:base 2> /dev/null)" != "" ]]); then
        echo "Docker container for PNI is built and tagged as pni:base"
    elif [ "$(docker images -q pni:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at pni:base"
    else
        echo "Oops! Unable to build Docker container for PNI"
    fi
    popd

    pushd $BASEDIR/Algorithms/PPCOR/
    docker build -t ppcor:base .
    if ([ $? = 0 ] && [[ "$(docker images -q ppcor:base 2> /dev/null)" != "" ]]); then
        echo "Docker container for PPCOR is built and tagged as ppcor:base"
    elif [ "$(docker images -q ppcor:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at ppcor:base"
    else
        echo "Oops! Unable to build Docker container for PPCOR"
    fi
    popd

    pushd $BASEDIR/Algorithms/SINGE/
    docker build -t singe:base .
    if ([ $? = 0 ] && [[ "$(docker images -q singe:base 2> /dev/null)" != "" ]]); then
        echo "Docker container for SINGE is built and tagged as singe:base"
    elif [ "$(docker images -q singe:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at singe:base"
    else
        echo "Oops! Unable to build Docker container for SINGE"
    fi
    popd

    pushd $BASEDIR/Algorithms/SCNS/
    docker build -t scns:base .
    if ([ $? = 0 ] && [[ "$(docker images -q scns:base 2> /dev/null)" != "" ]]); then
        echo "Docker container for SCNS is built and tagged as scns:base"
    elif [ "$(docker images -q scns:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at scns:base"
    else
        echo "Oops! Unable to build Docker container for SCNS"
    fi
    popd

    pushd $BASEDIR/Algorithms/SCODE/
    docker build -t scode:base .
    if ([ $? = 0 ] && [[ "$(docker images -q scode:base 2> /dev/null)" != "" ]]); then
        echo "Docker container for SCODE is built and tagged as scode:base"
    elif [ "$(docker images -q scode:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at scode:base"
    else
        echo "Oops! Unable to build Docker container for SCODE"
    fi
    popd

    pushd $BASEDIR/Algorithms/SCRIBE/
    docker build -t scribe:base .
    if ([ $? = 0 ] && [[ "$(docker images -q scribe:base 2> /dev/null)" != "" ]]); then
        echo "Docker container for SCRIBE is built and tagged as scribe:base"
    elif [ "$(docker images -q scribe:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at scribe:base"
    else
        echo "Oops! Unable to build Docker container for SCRIBE"
    fi
    popd

    pushd $BASEDIR/Algorithms/SINCERITIES/
    docker build -t sincerities:base .
    if ([ $? = 0 ] && [ "$(docker images -q sincerities:base 2> /dev/null)" != "" ]); then
        echo "Docker container for SINCERITIES is built and tagged as sincerities:base"
    elif [ "$(docker images -q sincerities:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at sincerities:base"
    else
        echo "Oops! Unable to build Docker container for SINCERITIES"
    fi
    popd

    pushd $BASEDIR/Algorithms/SCSGL/
    docker build -t scsgl:base .
    if ([ $? = 0 ] && [[ "$(docker images -q scsgl:base 2> /dev/null)" != "" ]]); then
        echo "Docker container for SCSGL is built and tagged as scsgl:base"
    elif [ "$(docker images -q scsgl:base 2> /dev/null)" != "" ]; then
        echo "Docker container failed to build, but an existing image exists at scsgl:base"
    else
        echo "Oops! Unable to build Docker container for SCSGL"
    fi
    popd
else
    echo "Pulling docker images from https://hub.docker.com/u/grnbeeline..."
    docker image pull grnbeeline/arboreto:base $VERBOSE_VALUE
    docker image pull grnbeeline/grisli:base $VERBOSE_VALUE
    docker image pull grnbeeline/grnvbem:base $VERBOSE_VALUE
    docker image pull grnbeeline/leap:base $VERBOSE_VALUE
    docker image pull grnbeeline/pidc:base $VERBOSE_VALUE
    docker image pull grnbeeline/ppcor:base $VERBOSE_VALUE
    docker image pull grnbeeline/scinge:base $VERBOSE_VALUE
    docker image pull grnbeeline/scns:base $VERBOSE_VALUE
    docker image pull grnbeeline/scode:base $VERBOSE_VALUE
    docker image pull grnbeeline/scribe:base $VERBOSE_VALUE
    docker image pull grnbeeline/sincerities:base $VERBOSE_VALUE
    docker image pull grnbeeline/singe:0.4.1 $VERBOSE_VALUE
fi
