#!/usr/bin/env python
# coding: utf-8

import os
import sys
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
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as line
import matplotlib.patches as patches

# local imports
import BLEval as ev

def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Generate plots from evaluation results.')

    parser.add_argument('-c','--config', default='config.yaml',
        help="Comma delimited configuration file(s) containing list of datasets "
              "algorithms and output specifications.\n")

    parser.add_argument('-a', '--auprc', action="store_true", default=False,
        help="Generate box plot of AUPRC values for evaluated algorithms.\n")
    
    parser.add_argument('-r', '--auroc', action="store_true", default=False,
        help="Generate box plot of AUROC values for evaluated algorithms.\n")
    
    parser.add_argument('-e', '--epr', action="store_true", default=False,
        help="Generate box plot of early precision values for evaluated algorithms.\n")

    parser.add_argument('-o', '--overview', action="store_true", default=False,
      help="Generate plot of AUPRC and early precision ratios relative to a random predictor.\n")

    return parser


def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts

def boxplot(evalConfigs, datasets, randValue, resTypeFile, resTypeName):
    
    # Set plot variables
    plt.rcParams.update({'font.size': 14})
    
    # Initialize plot window
    f, axes = plt.subplots(len(datasets), 1, figsize=(10, 7*len(datasets)))
    
    for i, dataset in enumerate(datasets):
        evalConfig = evalConfigs[i]
        
        # Read output file containing AUROC values
        DF = pd.read_csv(str(evalConfig.output_settings.base_dir) + '/' \
                           + str(evalConfig.input_settings.datadir).split("inputs")[1] + '/' \
                           + str(evalConfig.output_settings.output_prefix) \
                           + '-' + resTypeFile + '.csv', header = 0, index_col = 0)
        
        DF = DF.T
        modifiedDF = pd.melt(DF, id_vars=DF.index, value_vars=DF.columns)
              
        if len(datasets) > 1:
            subax = axes[i]
        else:
            subax= axes
        
        # Indicate random predictors value
        ax = sns.lineplot(y = randValue, 
                          x=range(-1,len(DF.columns)+1), ax = subax)
        ax.lines[0].set_linestyle("--")
        ax.lines[0].set_color("gray")
        
        # Plot the AUROC values as a box plot
        ax = sns.boxplot(y='value',x ='variable', data=modifiedDF, fliersize = 0,
                        palette=sns.color_palette("Set1"),
                        ax = subax)
        sns.swarmplot(y='value',x ='variable', data=modifiedDF,
                         alpha = 0.5, palette=['k'],
                         dodge = True,
                         ax = subax)
        
        ax.set_ylim([0.0,1])
        ax.set_ylabel(resTypeName, fontsize = 18)
        
        ax.set_title(resTypeName +' values for ' + dataset +' Network', fontsize = 18)
        plt.tight_layout()
    
    subax.set_xlabel('Algorithm', fontsize = 18)
               
    plt.savefig(str(evalConfig.output_settings.base_dir) + '/' \
                 + str(evalConfig.input_settings.datadir).split("inputs")[1] + '/' \
                 + str(evalConfig.output_settings.output_prefix) \
                 + '-boxplot-' + resTypeFile + '.pdf', dpi = 300)

def COplot(inputDF, width = 12, height = 7, randValues = [], refValues = [], rangeValues = [], shape = [], 
            palettes = [], text = [], levels = [], rotation = [], dispValues = [], switch = [], txtRange = (0,1)):

    plt.rcParams.update({'font.size': 12})

    levls = levels
    rowNames = inputDF.index
    maxRows = len(inputDF.index)
    maxCols = len(inputDF.columns)
    pad = 2
    fSize = (width,height)

    f = plt.figure(figsize=fSize)

    
    ax = plt.gca()
    ax.set_yticks(np.arange(0,maxRows+pad))
    ax.set_xticks(np.arange(0,maxCols+pad))
    

    Indices = [""] + list(rowNames)  +  [""]
    ax.set_yticklabels(Indices)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    alt = False
    for rowIdx in range(len(rowNames)):
        if alt:
            ptch = patches.Rectangle((0,rowIdx+0.5),
                               width = maxCols + 1,
                               height = 1,  
                               edgecolor = (1,1,1),
                               facecolor = (0.9,0.9,0.9),)

        else:
            ptch = patches.Rectangle((0,rowIdx+0.5),
                   width = maxCols,
                   height = 1,  
                   edgecolor = (1,1,1),
                   facecolor = (1,1,1),)
        alt = not(alt)
        ax.add_artist(ptch)

    colCnt = 0
    for levlIdx in range(len(levls)):
        randValue = randValues[levlIdx]
        refValue = refValues[levlIdx]
        ranges = rangeValues[levlIdx]

        colStart = colCnt
        plt.text(colStart+2, maxRows + 2, 
                        levls[levlIdx], fontsize=18, rotation=0,
                        ha="center", va="center",
                        bbox=dict(boxstyle="round",
                        ec=(1,1,1), fc=(1,1,1)))
        
        colNames = inputDF[levels[levlIdx]].columns

        for colIdx in range(len(colNames)):
                colCnt += 1

                plt.text(colStart + colIdx + 1, maxRows + 1, colNames[colIdx],
                         fontsize=18, rotation=rotation[levlIdx],
                         ha="center", va="center",
                         bbox=dict(boxstyle="round",
                         ec=(1,1,1,0),
                         fc=(1,1,1,0)))

                for rowIdx in range(len(rowNames)):

                    value = inputDF.loc[rowNames[rowIdx],levls[levlIdx]][colNames[colIdx]]
                    if text[levlIdx]:
                        if shape[levlIdx] == 'arrow':                    

                            valList = ['\u2197','\u2198','\u003d']

                            if value > 1:
                                txt = valList[0]
                                col = '#4CAF50'
                                fSize = 24

                            elif value < 1:
                                txt = valList[1]
                                col = '#E53935'
                                fSize = 24

                            elif value == 1:
                                txt = valList[2]
                                col = '#2E4053'
                                fSize = 24

                            else:
                                # Missing values?
                                txt = '-'
                                col = '#2E4053'
                                fSize = 24

                            plt.text(colStart+colIdx+1, rowIdx+1, 
                                txt, fontsize= fSize, rotation=0,
                                ha="center", va="center", color = col,
                                bbox=dict(boxstyle="round", 
                                ec=(1,1,1,0), fc=(1,1,1,0)))

                        else:
                            fSize = 23

                            if value < 1:
                                txt = round(value,2)
                            else:
                                txt = int(value)

                            col = '#2E4053'
                            fSize = 11
                            plt.text(colStart+colIdx+1, rowIdx+1, 
                                txt, fontsize= fSize, rotation=0,
                                ha="center", va="center", color = col,
                                bbox=dict(boxstyle="round", 
                                ec=(1,1,1,0), fc=(1,1,1,0)))              
                    
                    if shape[levlIdx] not in ['text','arrow']:
                        Oldvalue = value
                        value1 = value
                        
                        if np.isnan(value):
                            txt = '-'
                            col = '#2E4053'
                            fSize = 24

                            plt.text(colStart+colIdx+1, rowIdx+1, 
                            txt, fontsize= fSize, rotation=0,
                            ha="center", va="center", color = col,
                            bbox=dict(boxstyle="round", 
                            ec=(1,1,1,0), fc=(1,1,1,0)))
                            continue
                        elif value < randValue:
                            col = '#2E4053'
                        elif value >= ranges[1]:
                            if ranges[1] > 0:
                                value = round(value,1)
                                rangesD = (round(inputDF[levls[levlIdx]][colNames[colIdx]].min(),1), 
                                          round(inputDF[levls[levlIdx]][colNames[colIdx]].max(),1))

                                value1 = (value-rangesD[0])/(rangesD[1]-rangesD[0])
                                col = palettes[levlIdx][int(np.floor((value1)*10))]
                                value = 1
                            else:
                                value = round(value,3)
                                
                                rangesD = (round(inputDF[levls[levlIdx]][colNames[colIdx]].min(),3), 
                                          round(inputDF[levls[levlIdx]][colNames[colIdx]].max(),3))

                                value1 = (value-rangesD[0])/(rangesD[1]-rangesD[0])
                                col = palettes[levlIdx][int(np.floor((value1)*10))]
                                value = 1
                        else:
                            col = palettes[levlIdx][int(np.floor((value)*10))]


                        if shape[levlIdx] == 'c':
                            circle1=patches.Circle((colStart+colIdx+1,rowIdx+1),
                                               radius = np.sqrt(value)/2.5,
                                               facecolor=col,
                                               edgecolor = 'k',)

                        elif shape[levlIdx] == 's':
                            value = value*0.8
                            circle1=patches.Rectangle((colStart+colIdx+1-(value/2),rowIdx+1-(value)/2),
                                       width = value,
                                       height = value,    
                                       facecolor=col,
                                       edgecolor = 'k',)
                            
                        
                        elif shape[levlIdx] == 'hm':
                            value = 1
                            circle1=patches.Rectangle((colStart+colIdx+1-(value/2),rowIdx+1-(value)/2),
                                       width = value-0.1,
                                       height = value,    
                                       facecolor = col,
                                       edgecolor = col,)
                            

                        elif shape[levlIdx] == 'rs':
                            value = value*0.8
                            if value <= 0.15:
                                value = 0.15
                            boxPad = 0.075
                            newVal = value - boxPad*2
                            circle1=patches.FancyBboxPatch((colStart+colIdx+1-(newVal/2),rowIdx+1-(newVal)/2),
                                       width = newVal,
                                       height = newVal,    
                                       facecolor=col,
                                       edgecolor = 'k',
                                       boxstyle=patches.BoxStyle("Round", pad=boxPad))
                                
                        elif shape[levlIdx] == 'w':
                            circle1=patches.Wedge((colStart+colIdx+1,rowIdx+1),
                                       r = 0.4,
                                       theta1 = 0,
                                       theta2 = round(value*360,2),    
                                       facecolor = col,
                                       edgecolor = 'k',)
                        elif shape[levlIdx] == 'b':
                            circle1=patches.Rectangle((colStart+colIdx+0.6,rowIdx+0.65),
                                       width = 0.75,
                                       height = value,    
                                       facecolor=col,
                                       edgecolor = 'k',)
                            
                        elif shape[levlIdx] == 'f':
                            # flat color
                            circle1=patches.Rectangle((colStart+colIdx+1,rowIdx+0.6),
                                       width = 1,
                                       height = 1,    
                                       facecolor = col,
                                       edgecolor = col,)


                        ax.add_artist(circle1)
                        
                        textCol = ['black','white']
                        if switch[levlIdx]:
                            textCol[0] = 'white'
                            textCol[1] = 'black'

                        if value1 >= txtRange[1]:
                            plt.text(colStart+colIdx+1, rowIdx+1, 
                            round(Oldvalue,1), fontsize= 13, rotation=0,
                            ha="center", va="center", color = textCol[0],
                            bbox=dict(boxstyle="round", 
                            ec=(1,1,1,0), fc=(1,1,1,0)))
                        if value1 <= txtRange[0]:
                            plt.text(colStart+colIdx+1, rowIdx+1, 
                            round(Oldvalue,1), fontsize= 13, rotation=0,
                            ha="center", va="center", color = textCol[1],
                            bbox=dict(boxstyle="round", 
                            ec=(1,1,1,0), fc=(1,1,1,0)))

    ax.yaxis.set_ticks_position('none') 
    ax.xaxis.set_ticks_position('none') 
    ax.set_xticklabels([])

def main():
    opts = parse_arguments()
    config_files = opts.config.split(',')

    datasets = []
    evalConfigs = []

    for config_file in config_files:
        with open(config_file, 'r') as conf:
            evalConfig = ev.ConfigParser.parse(conf)
            evalConfigs.append(evalConfig)
            datasets.append(str(os.path.basename(evalConfig.input_settings.datadir)))

    randPredictor = {}
    # To compute network density for a dataset
    for i, dataset in enumerate(datasets):
        evalConfig = evalConfigs[i]
        
        # Read the reference network
        trueEdgesDF = pd.read_csv(str(evalConfig.input_settings.datadir) + '/' + evalConfig.input_settings.datasets[0]['name'] \
                               + '/' + evalConfig.input_settings.datasets[0]['trueEdges'],
                               header = 0, index_col = None)
        # Remove self-edges in network density computation, if any
        trueEdgesDF = trueEdgesDF[trueEdgesDF.Gene1 != trueEdgesDF.Gene2]
        # remove duplicated edges, if any
        trueEdgesDF.drop_duplicates(inplace = True)
        
        # Compute number of possible directed edges
        numGenes = len(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]))
        numPossibleEdges = numGenes*(numGenes-1)
        randPred = trueEdgesDF.shape[0]/numPossibleEdges
        print("Early Precision of a random predictor for ", dataset ," network (excluding self loops) is: %.2f" %(randPred))
        randPredictor[dataset] = randPred

    # Generate AUPRC boxplot     
    if (opts.auprc):
        print('\n\nGenerating AUPRC boxplot...')
        boxplot(evalConfigs, datasets, randPredictor[dataset], 'AUPRC', 'AUPRC')

    # Generate AUPRC boxplot     
    if (opts.auroc):
        print('\n\nGenerating AUROC boxplot...')
        boxplot(evalConfigs, datasets, 0.5, 'AUROC', 'AUROC')
        
    # Generate AUPRC boxplot     
    if (opts.epr):
        print('\n\nGenerating Early Precision boxplot...')
        boxplot(evalConfigs, datasets, randPredictor[dataset], 'EPr', 'Early Precision')

    # Generate overview plot     
    if (opts.overview):
        print('\n\nGenerating overview plot...')

        #TODO add option to select overview columns
        #resType = ['AUPRC Ratio','Early Precision Ratio']
        #resTypeFileName = [ 'AUPRC', 'EPr']
        resType = ['AUPRC Ratio', 'Stability Across Datasets']
        resTypeFileName = [ 'AUPRC', 'Jaccard']
        
        # Store results from each dataset in a dictionary and
        # then convert it to a dataframe.
        Res = {}
        
        # Initialize dataframe to store results
        multIndTuple = []
        for res in resType:
            for dataset in datasets:
                multIndTuple.append((res,dataset))
                
        ResDF = pd.read_csv(str(evalConfigs[0].output_settings.base_dir) + '/' \
                              + str(evalConfigs[0].input_settings.datadir).split("inputs")[1] + '/' \
                              + str(evalConfigs[0].output_settings.output_prefix) \
                              + '-' + resTypeFileName[0] + '.csv', header = 0, index_col = 0)
        algs = ResDF.index

        overviewDF = pd.DataFrame(index = algs, columns = pd.MultiIndex.from_tuples(multIndTuple))
        
        for i, res in enumerate(resType):
            for j, dataset in enumerate(datasets):
                evalConfig = evalConfigs[j]
                
                ResDF = pd.read_csv(str(evalConfig.output_settings.base_dir) + '/' \
                                    + str(evalConfig.input_settings.datadir).split("inputs")[1] + '/' \
                                    + str(evalConfig.output_settings.output_prefix) \
                                    + '-' + resTypeFileName[i] + '.csv', header = 0, index_col = 0)
                
                ResDF = ResDF.T
                for alg in algs:
                    # If early precision, compute the Early Precision Ratio (EPR)
                    # by dividing the values by that of a random predictor
                    if res == 'Early Precision Ratio' or res == 'AUPRC Ratio':
                        overviewDF.loc[alg][res,dataset] =  ResDF[alg].median()/randPredictor[dataset]
                    elif res == 'Stability Across Datasets' or res == 'Spearman':
                        overviewDF.loc[alg][res,dataset] =  ResDF.iloc[0][alg]
                    else:
        
                        overviewDF.loc[alg][res,dataset] =  ResDF[alg].median()
        
        overviewDF = overviewDF.loc[overviewDF['AUPRC Ratio'].median(axis='columns').sort_values().index]
        
        pale1 = sns.cubehelix_palette(11, dark=1, light = 0)
        pale2 = sns.color_palette("plasma_r",11)
        pale3 = sns.color_palette("viridis",11)
        pale3h = sns.cubehelix_palette(11, reverse = True)#[2:]
        pale3r = sns.color_palette("magma",11)
        pale4 = sns.color_palette("magma_r",13)[:-2]
        plt.rcParams["font.size"] = 14
        pale5 = sns.cubehelix_palette(rot=-0.3,n_colors = 11)
        pale5 = sns.color_palette('cividis_r', n_colors = 12)[:-1]
        COplot(overviewDF, width = 5, height = 7, randValues = [1, 0, 0, 0], 
                refValues = [1,0,0,0], rangeValues = [(1,1.001),(1.00,1.00),(0,0),(1,3)], 
                shape = ['hm','hm','hm','arrow'], palettes = [pale3, pale3h, pale5 ,pale3r],
                text = [False, False, False, True], 
                levels = resType,
                rotation = [0,0,0,0], dispValues =  [True, False, True, False], 
                switch = [False, False, True, False], txtRange = (0,1))
        
        plt.tight_layout()
        plt.savefig(str(evalConfig.output_settings.base_dir) + '/' \
                        + str(evalConfig.input_settings.datadir).split("inputs")[1] + '/' \
                        + str(evalConfig.output_settings.output_prefix) + '-overview.pdf')

if __name__ == '__main__':
  main()
