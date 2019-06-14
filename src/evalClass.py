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
import src.plotCurves as pc
import src.signedEPR as sEPR
import src.analyzeTopk as tk
import src.parseTime as pTime
import src.computeMotifs as cm
import src.computeEarlyPR as EPR
import src.computeSpearman as cSpear

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



class Evaluation(object):
    '''
    The Evaluation object is created by parsing a user-provided configuration
    file. Its methods provide for further processing its inputs into
    a series of jobs to be run, as well as running these jobs.
    '''

    def __init__(self,
            input_settings: InputSettings,
            output_settings: OutputSettings) -> None:

        self.input_settings = input_settings
        self.output_settings = output_settings

        
    def findSpearman(self):

        '''
        Finds the Spearman's correlation coefficient between consecutive simulations
        of the same algorithm, with the same initial parameters
        Saves the coefficients to a csv file, to be read into a Pandas dataframe
        :return:
        None
        '''
        corrDF = {}
        corrDF['Median'] = {}
        corrDF['MAD'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        print(outDir)
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                corrDF['Median'][algo[0]],corrDF['MAD'][algo[0]] = cSpear.Spearman(self, algo[0])
            
        return pd.DataFrame(corrDF)

    def findJaccard(self):

        '''
        Computes the Jaccard Index.
        Saves the coefficients to a csv file, to be read into a Pandas dataframe.
        :return:
        None
        '''
        JaccDF = {}
        JaccDF['Median'] = {}
        JaccDF['MAD'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        print(outDir)
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                JaccDF['Median'][algo[0]], JaccDF['MAD'][algo[0]] = self.computeJaccard(algo[0])
            
        return pd.DataFrame(JaccDF)
                 
    def findEarlyPr(self):

        '''
        Computes the Early Precision values.
        Saves the coefficients to a csv file, to be read into a Pandas dataframe.
        :return:
        None
        '''
        earlyPR = {}
        earlyPR['Early Precision'] = {}
        Eprec = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        print(outDir)
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                earlyPR['Early Precision'][algo[0]], Eprec[algo[0]] = EPR.EarlyPR(self, algo[0])
        pd.DataFrame(Eprec).to_csv(outDir+'EarlyPRFull.csv')
        return pd.DataFrame(earlyPR)
    
    def findSignedPr(self):

        '''
        Computes the Early Precision values for activation and inhibitory edges.
        Saves the coefficients to a csv file, to be read into a Pandas dataframe.
        :return:
        None
        '''
        sEPRDict = {}
        sEPRDict['Activation Precision'] = {}
        sEPRDict['Inhibition Precision'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                sEPRDict['Activation Precision'][algo[0]], sEPRDict['Inhibition Precision'][algo[0]] = sEPR.signedEPR(self, algo[0])
            
        return pd.DataFrame(sEPRDict)
    
         

    def findAUC(self):

        '''
        Plot PR and ROC curves for each dataset
        for all the algorithms
        '''
        AUPRCDict = {}
        AUROCDict = {}
        uAUPRCDict = {}
        uAUROCDict = {}
        
        for dataset in tqdm(self.input_settings.datasets, 
                            total = len(self.input_settings.datasets), unit = " Dataset"):
            
            AUPRC, AUROC = pc.PRROC(dataset, self.input_settings, directed = True, selfEdges = False)
            uAUPRC, uAUROC = pc.PRROC(dataset, self.input_settings, directed = False, selfEdges = False)
            
            AUPRCDict[dataset['name']] = AUPRC
            AUROCDict[dataset['name']] = AUROC
            uAUPRCDict[dataset['name']] = uAUPRC
            uAUROCDict[dataset['name']] = uAUROC
            
        return AUPRCDict, AUROCDict, uAUPRCDict, uAUROCDict
    
    def findMotifs(self):

        '''
        compute motifs
        '''
        
        for dataset in tqdm(self.input_settings.datasets, 
                            total = len(self.input_settings.datasets), unit = " Dataset"):
            cm.motifAnalysis(dataset, self.input_settings)
    
    def analyzeTopK(self):
        '''
        compute motifs
        '''
        
        for dataset in tqdm(self.input_settings.datasets, 
                            total = len(self.input_settings.datasets), unit = " Dataset"):
            tk.outputAnalysis(dataset, self.input_settings)
    
    def parseTimes(self):
        """
        Parse time output for each
        algorithm-dataset combination.
        
        Returns:
        -------
        TimeDict: 
            A dictionary of times for all dataset-algorithm combinations
        """
        TimeDict = dict()

        for dataset in self.input_settings.datasets:
            timevals  = pTime.getTime(self, dataset)
            TimeDict[dataset["name"]] = timevals

        return TimeDict
    
class ConfigParser(object):
    '''
    Define static methods for parsing a config file that sets a large number
    of parameters for the pipeline.
    '''
    @staticmethod
    def parse(config_file_handle) -> Evaluation:
        config_map = yaml.load(config_file_handle)
        return Evaluation(
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