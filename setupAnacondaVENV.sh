#!/bin/bash

# echo "If this is your first time installing Anaconda, you will need to run 'conda init bash' before running this script."
# echo "Run this command and then *restart the shell* for changes to take effect. Then run this script. Thank you!"
# read -p "Press Enter to continue..."

# Set-up Anaconda virtual environment
echo "Setting up Anaconda Python virtual environment..."

# You must have a root Python environment installed in Conda equivalent to 3.7.1 or higher because of an Anaconda bug.
# See https://github.com/ContinuumIO/anaconda-issues/issues/6698 as an example issue
conda install python # ^ Reasoning for this line
conda env create -f environment.yml 
eval "$(conda shell.bash hook)"
conda activate BEELINE

# Install the PRROC package for computing area under PR curve
# TODO: Write the PRROC AUC function and make it BEELINE package without using rpy2

R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/PRROC/PRROC_1.3.1.tar.gz', type = 'source')"
