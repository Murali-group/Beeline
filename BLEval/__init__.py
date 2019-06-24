import os
import yaml
import argparse
import itertools
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import multiprocessing
from pathlib import Path
import concurrent.futures
from itertools import permutations
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from networkx.convert_matrix import from_pandas_adjacency


# local imports
from BLEval.computeAUC import PRROC
from BLEval.parseTime import getTime
from BLEval.computeJaccard import Jaccard
from BLEval.computeSpearman import Spearman
from BLEval.computeNetMotifs import Motifs
from BLEval.computeEarlyPrec import EarlyPrec
from BLEval.computeSignedEPrec import signedEPrec


class InputSettings(object):
    '''
    Structure for storing the names of input files.
    '''

    def __init__(self,
            datadir, datasets, algorithms) -> None:

        self.datadir = datadir
        self.datasets = datasets
        self.algorithms = algorithms


class OutputSettings(object):
    '''
    Structure for storing the names of directories that output should
    be written to.
    '''

    def __init__(self, base_dir, output_prefix: Path) -> None:
        self.base_dir = base_dir
        self.output_prefix = output_prefix




class BLEval(object):
    '''
    The BEELINE Evaluation object is created by parsing a user-provided configuration
    file. Its methods provide for further processing its inputs into
    a series of jobs to be run, as well as running these jobs.
    '''

    def __init__(self,
            input_settings: InputSettings,
            output_settings: OutputSettings) -> None:

        self.input_settings = input_settings
        self.output_settings = output_settings


    def computeAUC(self, directedFlag = True):

        '''
        Plot PR and ROC curves for each dataset
        for all the algorithms
      
        Returns:
        -------
        '''
        AUPRCDict = {}
        AUROCDict = {}

        for dataset in tqdm(self.input_settings.datasets, 
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            
            AUPRC, AUROC = PRROC(dataset, self.input_settings, 
                                    directed = directedFlag, selfEdges = False, plotFlag = False)
            AUPRCDict[dataset['name']] = AUPRC
            AUROCDict[dataset['name']] = AUROC
            
        return AUPRCDict, AUROCDict
    

    def parseTime(self):
        """
        Parse time output for each
        algorithm-dataset combination.
        
        Returns:
        -------
        TimeDict: 
            A dictionary of times for all dataset-algorithm combinations
        """
        TimeDict = dict()

        for dataset in tqdm(self.input_settings.datasets, 
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            timevals  = getTime(self, dataset)
            TimeDict[dataset["name"]] = timevals

        return TimeDict

    def computeJaccard(self):

        '''
        Computes the Jaccard Index.
        Saves the coefficients to a csv file, to be read into a Pandas dataframe.
      
        Returns:
        -------
        '''
        JaccDF = {}
        JaccDF['Jaccard Median'] = {}
        JaccDF['Jaccard MAD'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                JaccDF['Jaccard Median'][algo[0]], JaccDF['Jaccard MAD'][algo[0]] = Jaccard(self, algo[0])
            
        return pd.DataFrame(JaccDF)     

    
    def computeSpearman(self):

        '''
        Finds the Spearman's correlation coefficient between consecutive simulations
        of the same algorithm, with the same initial parameters
        Saves the coefficients to a csv file, to be read into a Pandas dataframe
        
        Returns:
        -------
        corrDF: A pandas dataframe containing the median and median absolute
        deviation of the correlation values.
        '''
        corrDF = {}
        corrDF['Spearman Median'] = {}
        corrDF['Spearman MAD'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                corrDF['Spearman Median'][algo[0]],corrDF['Spearman MAD'][algo[0]] = Spearman(self, algo[0])
            
        return pd.DataFrame(corrDF)



    def computeNetMotifs(self):

        '''
        Finds network motifs such as FFL, FBL, and MI in the
        predicted network.
      
        Returns:
        -------
        '''
        FFLDict = {}
        FBLDict = {}
        MIDict = {}
        for dataset in tqdm(self.input_settings.datasets, 
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            FBLDict[dataset["name"]], FFLDict[dataset["name"]], MIDict[dataset["name"]] = Motifs(dataset, self.input_settings)
        
        return FBLDict, FFLDict, MIDict
                 
    def computeEarlyPrec(self):

        '''
        Computes the Early Precision values.
        Saves the coefficients to a csv file, to be read into a Pandas dataframe.
      
        Returns:
        -------
        '''
        earlyPR = {}
        earlyPR['Early Precision'] = {}
        Eprec = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                earlyPR['Early Precision'][algo[0]], Eprec[algo[0]] = EarlyPrec(self, algo[0])
        return pd.DataFrame(Eprec).T


    
    def computeSignedEPrec(self):

        '''
        Computes the Early Precision values for activation and inhibitory edges.
        Saves the coefficients to a csv file, to be read into a Pandas dataframe.
      
        Returns:
        -------
        '''
        sEPRDict = {}
        sEPRDict['Median EPR Activation'] = {}
        sEPRDict['Median EPR Inhibition'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                sEPRDict['Median EPR Activation'][algo[0]], sEPRDict['Median EPR Inhibition'][algo[0]] = signedEPrec(self, algo[0])
            
        return pd.DataFrame(sEPRDict)
    

class ConfigParser(object):
    '''
    Define static methods for parsing a config file that sets a large number
    of parameters for the pipeline.
    '''
    @staticmethod
    def parse(config_file_handle) -> BLEval:
        config_map = yaml.load(config_file_handle)
        return BLEval(
            ConfigParser.__parse_input_settings(
                config_map['input_settings']),
            ConfigParser.__parse_output_settings(
                config_map['output_settings']))

    @staticmethod
    def __parse_input_settings(input_settings_map) -> InputSettings:
        input_dir = input_settings_map['input_dir']
        dataset_dir = input_settings_map['dataset_dir']
        datasets = input_settings_map['datasets']

        return InputSettings(
                Path(input_dir, dataset_dir),
                datasets,
                ConfigParser.__parse_algorithms(
                input_settings_map['algorithms']))


    @staticmethod
    def __parse_algorithms(algorithms_list):
        algorithms = []
        for algorithm in algorithms_list:
                combos = [dict(zip(algorithm['params'], val))
                    for val in itertools.product(
                        *(algorithm['params'][param]
                            for param in algorithm['params']))]
                for combo in combos:
                    algorithms.append([algorithm['name'],combo])
            

        return algorithms

    @staticmethod
    def __parse_output_settings(output_settings_map) -> OutputSettings:
        output_dir = Path(output_settings_map['output_dir'])
        output_prefix = Path(output_settings_map['output_prefix'])

        return OutputSettings(output_dir,
                             output_prefix)
