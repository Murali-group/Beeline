"""
Sfaira Data Loader (:mod:`SfairaData`) module contains the following main class:


- :class:`SfairaData.SfairaData` and two additional classes used in the definition of SfairaData class
- :class:`SfairaData.SfairaSettings`
- :class:`SfairaData.ConfigParser`


"""


import yaml
import argparse
import itertools
from collections import defaultdict
from glob import glob
import pathlib
from pathlib import Path
import concurrent.futures
from typing import Dict, List
import multiprocessing
from multiprocessing import Pool, cpu_count
import concurrent.futures
import os
import pandas as pd
import sfaira
# import importlib_metadata
# import h5py
import scanpy as sc
import anndata
from andata import read_h5ad
import zipfile
import gzip
import tarfile
import urllib.request
import GEOparse




class SfairaSettings(object):
    '''    
    The class for storing the names of directories that datasets should
    be downloaded to and the features to filter subsets.
    This initilizes an SfairaSettings object based on the following parameters.
   
    :param base_dir: sfaira root directory, typically 'inputs/sfaira'
    :type base_dir: str
    :param filterss:   List of key-value pairs to filter subsets
    :type filters: list
    '''
    def __init__(self, base_dir, subsets) -> None:
        self.base_dir = base_dir
        self.subsets = subsets




class SfairaData(object):
    '''
    The SfairaData object is created by parsing a user-provided configuration
    file. Its methods provide for further processing its datasets into
    a series of jobs to be run, as well as running these jobs.
    '''
    def __init__(self,
                sfaira_settings: SfairaSettings) -> None:
        self.sfaira_settings = sfaira_settings




    def sfairaLoader(self):
        '''
        Download specified datasets from sfaira data repository.
   
        :returns:
            Folders containing scRNA-seq datasets.
        '''
        # Filter subset by looping all the key-value pairs
        for i in range(len(self.sfaira_settings)):
            basedir = self.sfaira_settings[i].base_dir
            datadir = os.path.join(basedir, 'raw')
            metadir = os.path.join(basedir, 'meta')
            cachedir = os.path.join(basedir, 'cache')
            ds = sfaira.data.Universe(data_path=datadir, meta_path=metadir, cache_path=cachedir)
            filters = self.sfaira_settings[i].subsets[i]


            ############# subset ############
            for key_val in filters:
                ds.subset(key=key_val, values=filters[key_val][:])


            ############# filtering ############
            # Download subsets to sepcifed folder and save in h5ad
            # Locate the cache folder
            cache_dir = os.path.join(dataset.path, "cache")


            # Move or copy the cache folder to the desired location
            new_cache_dir = "/path/to/your/desired/location/cache"
            shutil.copytree(cache_dir, new_cache_dir)  
 
            ############# download ############
            ds.download()
            ds.load()
            # Load the h5ad file from the cache folder
            # h5ad_path = os.path.join(new_cache_dir, f"{dataset.id}.h5ad")
            # adata = read_h5ad(h5ad_path)


   
    def csvConverter(self):
        '''
        Loop all files in multi-level subfolders
        '''
        for i in range(len(self.sfaira_settings)):
            subset_dir = self.sfaira_settings[i].base_dir
            SfairaData.__folderProcess(subset_dir)


    def __folderProcess(path):
        '''
        Convert files to csv files.
        '''
        if os.path.isdir(path):
            # The file is a folder, so loop through its contents
            for file in os.listdir(path):
                # Construct the full path to the item
                file_path = os.path.join(path, file)
               
                if os.path.isdir(file_path):
                    # The file is a folder, so recurse into it
                    SfairaData.__folderProcess(file_path)
                else:
                    # The file is not a folder, convert the file to csv file
                    # compressed files
                    if file.endswith('.zip'):
                        zipfile.ZipFile(file_path, "r").extractall(os.path.dirname(os.path.abspath(file_path)))
                    elif file.endswith('.gz'):
                        with gzip.open(file_path, 'rb') as f_in:
                            with open(file_path.replace("gz", "csv"), 'wb') as f_out:
                                f_out.write(f_in.read())
                    elif file.endswith('.tar'):
                        tarfile.open(file_path, "r").extractall(os.path.dirname(os.path.abspath(file_path)))


                    # h5ad files
                    if file.endswith('.h5ad'):
                        # Load h5ad file
                        adata = anndata.read_h5ad(file_path)
                        # Convert AnnData object to pandas DataFrame as gene x cell
                        df = pd.DataFrame(adata.X.todense(), index=adata.obs_names, columns=adata.var_names).T
                        # Write DataFrame to CSV file
                        file = file_path.replace("h5ad", "csv")
                        df.to_csv(file)
                                       
                        # # acc.gci files
                        # if file.startswith('acc.cgi'):
                        #     # Get the GEO accession number
                        #     accession_number = file.split('=')[1]
                        #
                        #     # Get data from python package GEOparse
                        #     # Download the metadata and expression data for the GEO accession number
                        #     gse = GEOparse.get_GEO(accession_number)
                        #     # Access the expression data for the dataset
                        #     expression_data = gse.table
                        #     # Save the expression data to a CSV file
                        #     file = accession_number + "csv"
                        #     expression_data.to_csv(file)
                        #
                        #     # Get data from URL
                        #     # Construct the URL for the GEO accession number
                        #     url = f'https://www.ncbi.nlm.nih.gov/geo/download/?acc={accession_number}&format=file'
                        #     # Set the file for the downloaded file
                        #     file = f'{accession_number}.RAW.tar'
                        #     # Download the file using the urllib module
                        #     urllib.request.urlretrieve(url, file)
 
class ConfigParser(object):
    '''
    The class define static methods for parsing and storing the contents
    of the config file.
    '''
    @staticmethod
    def parse(config_file_handle) -> SfairaData:
        '''
        A method for parsing the config .yaml file.
       
        :param config_file_handle: Name of the .yaml file to be parsed
        :type config_file_handle: str
       
        :returns:
            An object of class :class:`SfairaData.SfairaData`.


        '''
        config_map = yaml.load(config_file_handle, Loader=yaml.Loader)
        return SfairaData(
            ConfigParser.__parse_sfaira_settings(
                config_map['sfaira_settings']))
   


    @staticmethod
    def __parse_sfaira_settings(sfaira_settings_map) -> SfairaSettings:
        '''
        A method for parsing and initializing sfaira data object.
        '''
        # Obtain the data directory
        sfaira_dir = sfaira_settings_map['data_dir']


        # Obtain the subdata directory
        SfairaSettingsDir = {}
        key_order = ["year", "organism", "organ", "assay_sc"]
        order = 0


        for x in ConfigParser.__parse_sfaira_features(sfaira_settings_map['subsets']):
            subset_dir = ""
            for key in key_order:
                if key in x:
                    # Replace sapce in string with hypen  
                    val = str(x[key]).replace(", ", "+").replace(",", "+").replace(" ", "-").replace("[", "").replace("]", "").replace("'", "")
                    if len(subset_dir) == 0:
                        subset_dir = subset_dir + val
                    else:                                
                        subset_dir = subset_dir + "_" + val  
            # Set SfairaSettings for each specification of subset in config yaml
            # print(subset_dir)
            SfairaSettingsDir[order] = SfairaSettings(Path(sfaira_dir, subset_dir),
                                                           ConfigParser.__parse_sfaira_features(sfaira_settings_map['subsets']))
            order = order + 1
        return SfairaSettingsDir
   


    @staticmethod
    def __parse_sfaira_features(sfaira_list):
        '''
        A method for parsing parameters that determine the subsets to be downloaded.
        '''      
        # Initilalize the list of subset values
        subsets = []
        # print(sfaira_list)


        # Parse contents of sfaira_list
        for x in sfaira_list:
            key_values = x['filters']
            subsets.append(key_values)
        return subsets  
