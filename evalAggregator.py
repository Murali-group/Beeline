#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import eval as ev
import numpy as np
import pandas as pd
from tqdm import tqdm
import networkx as nx
import multiprocessing
from pathlib import Path
import concurrent.futures
from scipy.stats import spearmanr
from itertools import permutations
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from networkx.convert_matrix import from_pandas_adjacency

# local imports
import src.plotCurves as pc
import src.analyzeTopk as tk
import src.computeMotifs as cm
import src.parseTime as pTime
import src.computeSpearman as cSpear
import src.singedEPR as sEPR
import src.computeEarlyPR as EPR
 


class Evaluation(object):
    '''
    The Evaluation object is created by parsing a user-provided configuration
    file. Its methods provide for further processing its inputs into
    a series of jobs to be run, as well as running these jobs.
    '''

    def __init__(self,
            input_settings: ev.InputSettings,
            output_settings: ev.OutputSettings) -> None:

        self.input_settings = input_settings
        self.output_settings = output_settings

        
    def find_correlations(self):

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

    def find_jaccard(self):

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
        Computes the Jaccard Index.
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
        Computes the Jaccard Index.
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
    
         

    def find_curves(self):

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
    
    def find_motifs(self):

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
    
    def time_runners(self):
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


def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('--config', default='config.yaml',
        help='Configuration file')

    return parser

def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts

def main():
    opts = parse_arguments()
    config_file = opts.config

    evalConfig = None

    with open(config_file, 'r') as conf:
        evalConfig = ev.ConfigParser.parse(conf)
        
    print('Post-Run Evaluation Summary started')
    evalSummarizer = Evaluation(evalConfig.input_settings, evalConfig.output_settings)
    outDir = str(evalSummarizer.output_settings.base_dir) + \
            str(evalSummarizer.input_settings.datadir).split("inputs")[1] + "/"+\
            str(evalSummarizer.output_settings.output_prefix) + "-"
    
    # Compute time taken
    TimeDict = evalSummarizer.time_runners()
    pd.DataFrame(TimeDict).to_csv(outDir+'Timescores.csv')
    
    # Compute and plot ROC, PRC and report AUROC, AUPRC
    AUPRCDict, AUROCDict, uAUPRCDict, uAUROCDict = evalSummarizer.find_curves()

    pd.DataFrame(AUPRCDict).to_csv(outDir+'AUPRCscores.csv')
    pd.DataFrame(AUROCDict).to_csv(outDir+'AUROCscores.csv')
    pd.DataFrame(uAUPRCDict).to_csv(outDir+'uAUPRCscores.csv')
    pd.DataFrame(uAUROCDict).to_csv(outDir+'uAUROCscores.csv')
    
    # Compute early recision for signed edges
    sPrDF = evalSummarizer.findSignedPr()
    sPrDF.to_csv(outDir + "SingedPr.csv")
    
     # Compute early precision
    ePRDF = evalSummarizer.findEarlyPr()
    ePRDF.to_csv(outDir + "EarlyPR.csv")

    # Compute Jaccard index
    jaccDict = evalSummarizer.find_jaccard()
    jaccDict.to_csv(outDir + "jaccData.csv")

    # Compute correlations
    corrDict = evalSummarizer.find_correlations()
    corrDict.to_csv(outDir + "corrData.csv")


    print('Evaluation complete...')


if __name__ == '__main__':
  main()
