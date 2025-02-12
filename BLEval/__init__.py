
"""
BEELINE Evaluation (:mod:`BLEval`) module contains the following
:class:`BLEval.BLEval` and three additional classes used in the
definition of BLEval class 

- :class:`BLEval.ConfigParser` 
- :class:`BLEval.InputSettings` 
- :class:`BLEval.OutputSettings`


"""
import os
import yaml
import argparse
import itertools
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import multiprocessing
from glob import glob
import pathlib
from pathlib import Path
import concurrent.futures
from itertools import permutations
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from networkx.convert_matrix import from_pandas_adjacency


# local imports
from BLEval.parseTime import getTime
from BLEval.computeDGAUC import PRROC
from BLEval.computeBorda import Borda
from BLEval.computeJaccard import Jaccard
from BLEval.computeSpearman import Spearman
from BLEval.computeNetMotifs import Motifs
from BLEval.computeEarlyPrec import EarlyPrec
from BLEval.computePathStats import pathAnalysis
from BLEval.computeSignedEPrec import signedEPrec


class InputSettings(object):
    '''
    The class for storing the names of input files.
    This initilizes an InputSettings object based on the
    following three parameters.
    

    :param datadir:   input dataset root directory, typically 'inputs/'
    :type datadir: str

    :param datasets:   List of dataset names
    :type datasets: list
        
    :param algorithms:   List of algorithm names
    :type algorithms: list
    '''

    def __init__(self, datadir, datasets, algorithms, use_embeddings=False) -> None:
        self.datadir = Path(datadir)
        self.datasets = datasets
        self.algorithms = algorithms
        self.use_embeddings = use_embeddings

        if use_embeddings:
            self.adjust_paths_for_embeddings()
    
    def adjust_paths_for_embeddings(self):
        for dataset in self.datasets:
            processed_path = 'processed_ExpressionData'
            dataset['exprData'] = str(Path(processed_path) / 'EmbeddingsData.csv')
            dataset['refNetwork'] = str(Path(processed_path) / 'refNetwork.csv')
    
    def get_true_edges_path(self, dataset_name):
        dataset = next((d for d in self.datasets if d['name'] == dataset_name), None)
        if dataset is None:
            raise ValueError(f"Dataset with name {dataset_name} not found in datasets list.")
        return self.datadir / dataset_name / 'processed_ExpressionData' / dataset['trueEdges']
        
class OutputSettings(object):
    '''    
    The class for storing the names of directories that output should
    be written to. This initializes an OutputSettings object based on the
    following two parameters.
    
    :param base_dir: output root directory, typically 'outputs/'
    :type base_dir: str
    :param output_prefix: A prefix added to the final output files.
    :type str:
    '''

    def __init__(self, base_dir, output_prefix: Path, datasets, use_embeddings=False) -> None:
        self.base_dir = Path(base_dir)
        self.output_prefix = output_prefix
        self.datasets = datasets

        if use_embeddings:
            self.adjust_paths_for_embeddings()

    def adjust_paths_for_embeddings(self):
        self.base_dir = self.base_dir / 'processed_ExpressionData'

    def get_output_path(self, dataset_name, file_name):
        dataset = next((d for d in self.datasets if d['name'] == dataset_name), None)
        if dataset is None:
            raise ValueError(f"Dataset with name {dataset_name} not found in datasets list.")
        
        # Construct the output path based on the required structure
        output_path = self.base_dir / dataset_name / 'processed_ExpressionData' / f"{dataset_name}-{file_name}"
        return output_path

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


    def computeAUC(self, directed = True):
        '''
        Computes areas under the precision-recall (PR) and
        and ROC plots for each algorithm-dataset combination.
      
        Parameters
        ----------
        directedFlag: bool
            A flag to specifiy whether to treat predictions
            as directed edges (directed = True) or 
            undirected edges (directed = False).
        
        :returns:
            - AUPRC: A dataframe containing AUPRC values for each algorithm-dataset combination
            - AUROC: A dataframe containing AUROC values for each algorithm-dataset combination
        '''
        AUPRCDict = {}
        AUROCDict = {}

        for dataset in tqdm(self.input_settings.datasets, 
                            total=len(self.input_settings.datasets), unit=" Datasets"):
            if self.input_settings.use_embeddings == True:
                true_edges_path = self.input_settings.get_true_edges_path(dataset['name'])        

            AUPRC, AUROC = PRROC(dataset, self.input_settings, 
                                directed=directed, selfEdges=False, plotFlag=False)
            AUPRCDict[dataset['name']] = AUPRC
            AUROCDict[dataset['name']] = AUROC
              
        AUPRC = pd.DataFrame(AUPRCDict)
        AUROC = pd.DataFrame(AUROCDict)
        return AUPRC, AUROC
    

    def parseTime(self):
        """
        Parse time output for each
        algorithm-dataset combination.
        
        :returns:
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
        Computes Jaccard Index between top-k edge predictions
        of the same algorithm.
        
        :returns:
            A dataframe containing the median and median absolute
            deviation of the Jaccard Index values of each algorithm
            on the given set of datasets.
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
        Finds the Spearman's correlation coefficient between
        the ranked edges of the same algorithm on the given
        set of datasets.
        
        :returns:
            A dataframe containing the median and median absolute
            deviation of the Separman's correlation values of each algorithm.
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
        For each algorithm-dataset combination, this function computes the network motifs such as 
        Feedforward loops, Feedback loops and 
        Mutual interactions in the predicted top-k network. It returns the ratio of network motif counts
        compared to their respective values in the reference network.
      
        :returns:
            - FBL: A dataframe containing ratios of number of Feedback loops
            
            - FFL: A dataframe containing ratios of number of Feedforward loops
            
            - MI: A dataframe containing ratios of number of Mutual Interactions
        '''
        FFLDict = {}
        FBLDict = {}
        MIDict = {}
        for dataset in tqdm(self.input_settings.datasets, 
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            FBLDict[dataset["name"]], FFLDict[dataset["name"]], MIDict[dataset["name"]] = Motifs(dataset, self.input_settings)
            
        FBL = pd.DataFrame(FBLDict)
        FFL = pd.DataFrame(FFLDict)
        MI = pd.DataFrame(MIDict)
        
        return FBL, FFL, MI
    
    
    def computePaths(self):
        '''
        For each algorithm-dataset combination, this function computes path lengths
        through TP edges and FP edges, returns statistics on path lengths.

        :returns:
            - pathStats: A dataframe path lengths in predicted network
        '''
        for dataset in tqdm(self.input_settings.datasets, 
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            pathAnalysis(dataset, self.input_settings)

                 
    def computeEarlyPrec(self):

        '''
        For each algorithm-dataset combination,
        this function computes the Early Precision values of the 
        network formed using the predicted top-k edges.
        
        :returns: 
            A dataframe containing the early precision values 
            for each algorithm-dataset combination.

        '''
        Eprec = {}
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                Eprec[algo[0]] = EarlyPrec(self, self.input_settings, algo[0])
        return pd.DataFrame(Eprec).T
        
    def computeSignedEPrec(self):

        '''
        For each algorithm-dataset combination,
        this function computes the Early Precision values separately 
        for the activation and inhibitory edges.
        
        :returns: 
            - A dataframe containing early precision for activation edges
            - A dataframe containing early precision for inhibitory edges
        '''
        sEPRDict = {}
        sEPRDict['EPrec Activation'] = {}
        sEPRDict['EPrec Inhibition'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                sEPrecDF = signedEPrec(self, algo[0])
                sEPRDict['EPrec Activation'][algo[0]] = sEPrecDF['+']
                sEPRDict['EPrec Inhibition'][algo[0]] = sEPrecDF['-']
        return(pd.DataFrame(sEPRDict['EPrec Activation']).T, pd.DataFrame(sEPRDict['EPrec Inhibition']).T)
    
    def computeBorda(self, selectedAlgorithms=None, aggregationMethod="average"):

        '''
        Computes edge ranked list using the Borda method for each dataset.
        Parameters
        ----------
        selectedAlgorithms: [str]
          List of algorithm names used to run borda method on selected
          algorithms. If nothing is provided, the function runs borda on
          all available algorithms.
        aggregationMethod: str
          Method used to aggregate rank in borda method. Available options are
          {‘average’, ‘min’, ‘max’, ‘first’}, default ‘average’
        :returns:
            None
        '''
        feasibleAlgorithmOptions = [algorithmName for algorithmName, _ in self.input_settings.algorithms]
        feasibleaggregationMethodOptions = ['average', 'min', 'max', 'first']

        selectedAlgorithms = feasibleAlgorithmOptions if selectedAlgorithms is None else selectedAlgorithms

        for a in selectedAlgorithms:
            if a not in feasibleAlgorithmOptions:
                print("\nERROR: No data available on algorithm %s. Please choose an algorithm from the following options: %s" % (a, feasibleAlgorithmOptions))
                return

        if aggregationMethod not in feasibleaggregationMethodOptions:
                print("\nERROR: Please choose an aggregation method algorithm from following options: " % feasibleaggregationMethodOptions)
                return

        Borda(self, selectedAlgorithms, aggregationMethod)

class ConfigParser(object):
    '''
    The class define static methods for parsing and storing the contents 
    of the config file that sets a that sets a large number of parameters 
    used in the BLEval.
    '''
    @staticmethod
    def parse(config_file_handle, use_embeddings=False) -> BLEval:
        '''
        A method for parsing the input .yaml file.
        
        :param config_file_handle: Name of the .yaml file to be parsed
        :type config_file_handle: str
        
        :returns: 
            An object of class :class:`BLEval.BLEval`.

        '''
        config_map = yaml.load(config_file_handle, Loader=yaml.SafeLoader)
        input_settings = ConfigParser.__parse_input_settings(config_map['input_settings'], use_embeddings)
        output_settings = ConfigParser.__parse_output_settings(config_map['output_settings'], input_settings.datasets, use_embeddings)
        return BLEval(input_settings, output_settings)
    
    @staticmethod
    def __parse_input_settings(input_settings_map, use_embeddings=False) -> InputSettings:
        '''
        A method for parsing and initializing 
        InputSettings object.
        '''
        input_dir = input_settings_map['input_dir']  # e.g., 'inputs'
        dataset_dir = input_settings_map['dataset_dir']  # e.g., 'example'
        
        if 'datasets' in input_settings_map: 
            datasets_specified = True
        else:
            datasets_specified = False

        if datasets_specified is False: 
             subfolder_dir = glob(os.path.join(input_dir, dataset_dir, "*/"), recursive = True)
             datasets = []
             for x in subfolder_dir:
                datasets.append({"name": pathlib.Path(x).name, 
                                 "exprData": "ExpressionData.csv", 
                                 "cellData": "PseudoTime.csv", 
                                 "trueEdges": "refNetwork.csv"})
        # If datasets specified, run the corresponding datasets
        else:
            datasets = input_settings_map['datasets']
            if datasets is None:
                print("Please specify input datasets!")
        

        # datasets = input_settings_map['datasets']

        datadir = Path(input_dir) / dataset_dir  # e.g., 'inputs/example'
        return InputSettings(datadir, datasets, ConfigParser.__parse_algorithms(input_settings_map['algorithms']), use_embeddings)

    @staticmethod
    def __parse_algorithms(algorithms_list):
        '''
        A method for parsing the list of algorithms
        that are being evaluated, along with
        any parameters being passed.
        
        Note that these parameters may not be
        used in the current evaluation, but can 
        be used at a later point.
        '''
        
        # Initilalize the list of algorithms
        algorithms = []
        
        # Parse contents of algorithms_list
        for algorithm in algorithms_list:
                combos = [dict(zip(algorithm['params'], val))
                    for val in itertools.product(
                        *(algorithm['params'][param]
                            for param in algorithm['params']))]
                for combo in combos:
                    algorithms.append([algorithm['name'],combo])
            

        return algorithms

    @staticmethod
    def __parse_output_settings(output_settings_map,datasets,use_embeddings=False) -> OutputSettings:
        '''
        A method for parsing and initializing 
        Output object.
        '''
        output_dir = Path(output_settings_map['output_dir'])
        output_prefix = Path(output_settings_map['output_prefix'])
        return OutputSettings(output_dir, output_prefix, datasets, use_embeddings)
