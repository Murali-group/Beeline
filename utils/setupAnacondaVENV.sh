#!/bin/bash

# echo "If this is your first time installing Anaconda, you will need to run 'conda init bash' before running this script."
# echo "Run this command and then *restart the shell* for changes to take effect. Then run this script. Thank you!"
# read -p "Press Enter to continue..."

set -ex # abandon script on error and print out lines run

# Set-up Anaconda virtual environment
echo "Setting up Anaconda Python virtual environment..."

BASEDIR="$(dirname "$(readlink -f "$0")")"

# You must have a root Python environment installed in Conda equivalent to 3.7.1 or higher because of an Anaconda bug.
# See https://github.com/ContinuumIO/anaconda-issues/issues/6698 as an example issue
ENV_NAME="BEELINE"

if conda info --env | grep -q "${ENV_NAME}"; then
    echo "Conda environment '${ENV_NAME}' already exists. Checking for updates."
    conda env update --file=$BASEDIR/environment.yml -n "${ENV_NAME}"
else
    echo "Creating Conda environment '${ENV_NAME}'..."
    conda env create --file=$BASEDIR/environment.yml
    echo "Conda environment '${ENV_NAME}' created."
fi

# Set up per-algorithm environments for local (no-Docker) execution
for ALG_ENV in BEELINE-PPCOR BEELINE-LEAP BEELINE-GRNBOOST2; do
    # Derive yml filename: e.g. BEELINE-PPCOR -> environment-ppcor.yml
    ALG_YML="$BASEDIR/environment-$(echo ${ALG_ENV#BEELINE-} | tr '[:upper:]' '[:lower:]').yml"
    if conda info --env | grep -q "${ALG_ENV}"; then
        echo "Conda environment '${ALG_ENV}' already exists. Checking for updates."
        conda env update --file="$ALG_YML" -n "${ALG_ENV}"
    else
        echo "Creating Conda environment '${ALG_ENV}'..."
        conda env create --file="$ALG_YML"
        echo "Conda environment '${ALG_ENV}' created."
    fi
done

