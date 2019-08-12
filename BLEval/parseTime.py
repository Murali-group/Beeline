import os
import argparse
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

def getTime(evalObject, dataset):
    """
    Return time taken for each of the algorithms
    in the evalObject on the dataset specified.
    The output stored in time.txt is parsed to
    obtain the CPU time.

    :param evalObject:   An object of the class :class:`BLEval.BLEval`
    :type evalObject: BLEval
    :param dataset:   Dataset name for which the time output must be parsed for each algorithm.
    :type dataset: str
      
    
    :returns: 
        A Dictionary of time taken by each of the algorithms, i.e., key is the algorithm name, and
        value is the time taken (in sec.).

    """

    
    # Set the output directory for a given dataset 
    # where the time.txt files are stored for each algorithm.
    outDir = str(evalObject.output_settings.base_dir) + \
             str(evalObject.input_settings.datadir).split("inputs")[1] + "/" + dataset["name"] + "/"
    
    # Initialize algorithm:time-taken dictionary
    algo_dict = dict()
    
    # Obtain the list of algorithms for which the 
    # time.txt files must be read in the outDir
    algos = evalObject.input_settings.algorithms
    
    # Repeat for each algorithm
    for algo in algos:
        # Check if single time.txt file
        path = outDir+algo[0]+"/time.txt"
        if Path(path).exists():
            # If so, parse the time.txt file
            # to obtain CPU time taken to run the
            # GRN algorithm on the given dataset
            time = parse_time_files(path)
        else:
            # If the algorithm was run on each
            # trajectory in the dataset separately, there
            # will be multiple time files, namely, 
            # time0.txt, time1.txt, and so on. 

            # In order to report the total time taken,
            # we sum up the time taken for running the alorithm
            # on each trajectory.
            
            # First, obtain the number of trajectories
            # by reading the 'cellData' file for the given
            # dataset.
            
            # path to pseduotime file (or 'cellData' file)
            PTFile = str(evalObject.input_settings.datadir.joinpath(dataset['name']+'/'+dataset['cellData']))
         
            # Read in the PTFile to obtain the number of trajectories
            # Which is equal to the number of columns.           

            PTData = pd.read_csv(PTFile, header = 0, index_col = 0)
            
            # TODO: Check if column names only correspond to PseudoTime for 
            # each trajectory!
            colNames = PTData.columns
            
            # Read each time_.txt file individually
            # and report the total time taken.
            timeTemp = 0
            for idx in range(len(colNames)):
                path = outDir+algo[0]+"/time"+str(idx)+".txt"
                timeTemp += parse_time_files(path)

            # If for some reason a time file is missing,
            # the parse_time_files function returns a -1.
            # We will ignore those algorithms with missing
            # time values for now, as shown below.
            if timeTemp >= 0:
                time =  timeTemp
            else:
                time = -1
                
        # If time files do not exist, skip reporting that algorithm
        if time == -1:
            print("Skipping time computation for ", algo[0], "on dataset", dataset["name"], "\n")
            continue

        # If time files do exist, add it to the dictionary
        algo_dict[algo[0]] = time

    return algo_dict

def parse_time_files(path):
    """
    Return time taken for each of the algorithms
    in the evalObject on the dataset specified.
    The output stored in time.txt is parsed to
    obtain the CPU time.

    :param path: Path to the time.txt file, or timex.txt file where x corresponds to the trajectory ID for a given algorithm-dataset combination.
    :type path: str
      
    :returns: 
        A float value corresponding to the time taken.
         
    """
    try:
        with open(path, "r") as f:
            lines = f.readlines()
            line = lines[1]
            time_val = float(line.split()[-1])
            
    # If file is not found, return -1.
    except FileNotFoundError:
        print("Time output " +path+" file not found, setting time value to -1\n")
        time_val = -1
        
    # If file is present but the file is empty, return -1.
    except ValueError:
        print("Algorithm running failed, setting time value to -1\n")
        time_val = -1

    return time_val
