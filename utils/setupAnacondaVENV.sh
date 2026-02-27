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

eval "$(conda shell.bash hook)"
conda activate BEELINE

