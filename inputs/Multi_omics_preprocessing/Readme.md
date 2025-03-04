This directory is focussed on preprocessing of files for running Multi-omics based algorithms like MICA and DeepMAPS.
Most of these datasets are very large and take significant time and compute (at least 16GB RAM is advisable) to process. 
These preprocessing steps are mainly done based on the requirements of the MICA dataset and using standard bioinformatics tools, if better preprocessing is possible please do update the same. 

1. hTcell.py : this script takes in downloaded RNA and ATAC files and generates the regulators and ATAC processed files
It integrates the two modalities after normalization adn filtering. It uses spearman rank correlation between peaks and genes to assign correlation scores to the genes
Finally, using API calls, it assigns the relevant gene ames  for the ucsc id names of the nearest transcription start site. This helps link the genes to their regulatory genes.


2. mRetina : 
   - process_peaks.sh : Here the preprocessing is one on RDS data but similar datasets with h5 formats could also be used with minor alterations to the file read statements
     This downloads the gtf annotation files, the input dataset (peaks list file) . Then it links the R script to read RDS files and later the python processing file for combining and annotating the processed files
     This also used bedtools to find intersecting genes in the peak region and closest gene to the peak. This method helps link the peak given to the relevant gene (regulatory gene and target gene relationship)
     Finally, multiple intermediate files are removed. 

   - process_peaks.R : this script reads the downloaded files and makes the necessary bed and csv files for the peaks matrix and RNA matrix respectively. 

   - final_process.py : this script fetches gene names for ensemble id names, combines the nearest gene data and the annotated pek data based on location on chromosomes. And finally saves the processed files. 


    chmod +x process_peaks.sh 

    ./process_peaks.sh


3. Lymph-Node Lymphoma dataset : This script downloads necessary files related to the lymphoma dataset including the JASPAR reference files.


4. Other possible datasets : More datasets to be used for multiomics based approaches to build gene regulatory networks have been listed and deatlied in this recent paper https://www.biorxiv.org/content/10.1101/2024.02.01.578507v1.full