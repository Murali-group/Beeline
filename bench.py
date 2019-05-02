#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import yaml
import argparse
import itertools
from collections import defaultdict
from pathlib import Path
import multiprocessing
from multiprocessing import Pool, cpu_count
import concurrent.futures
from typing import Dict, List
from src.runner import Runner
import os, sys
#from src.plotCurves import EvalCurves
from src.eval import Eval
import pandas as pd
import fcntl

class InputSettings(object):
    def __init__(self,
            datadir, datasets, algorithms) -> None:
        
        self.datadir = datadir
        self.datasets = datasets
        self.algorithms = algorithms


class OutputSettings(object):
    '''
    Structure for storing the names of directories that output should
    be written to
    '''

    def __init__(self, base_dir: Path) -> None:
        self.base_dir = base_dir

        
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
        self.runners: Dict[int, Runner] = self.__create_runners()


    def __create_runners(self) -> Dict[int, List[Runner]]:
        '''
        Instantiate the set of runners based on parameters provided via the
        configuration file. Each runner is supplied an interactome, collection,
        the set of algorithms to be run, and graphspace credentials, in
        addition to the custom parameters each runner may or may not define.
        '''
        
        runners: Dict[int, Runner] = defaultdict(list)
        order = 0
        for dataset in self.input_settings.datasets:
            for runner in self.input_settings.algorithms:
                data = {}
                data['name'] = runner[0]
                data['params'] = runner[1]
                data['inputDir'] = Path.cwd().joinpath(self.input_settings.datadir.joinpath(dataset['name']))
                data['exprData'] = dataset['exprData']
                data['cellData'] = dataset['cellData']
                if 'should_run' in data['params'] and \
                        data['params']['should_run'] is False:
                    print("Skipping %s" % (data['name']))
                    continue
                #else:
                #    print("Running %s" % (data['name']))

                runners[order] = Runner(data)
                order += 1            
        return runners


    def execute_runners(self, parallel=False, num_threads=1):
        '''
        Run each of the algorithms
        '''

        base_output_dir = self.output_settings.base_dir

        batches =  self.runners.keys()

        for batch in batches:
            if parallel==True:
                executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
                futures = [executor.submit(runner.run, base_output_dir)
                    for runner in self.runners[batch]]
                
                # https://stackoverflow.com/questions/35711160/detect-failed-tasks-in-concurrent-futures
                # Re-raise exception if produced
                for future in concurrent.futures.as_completed(futures):
                    future.result()
                executor.shutdown(wait=True)
            else:
                for runner in self.runners[batch]:
                    runner.run(output_dir=base_output_dir)


    def evaluate_runners(self):
        '''
        Write AUPRC and AUROC summary statistics for all algorithms
        '''
        for dataset in self.input_settings.datasets:              
            #EvalCurves(dataset, self.input_settings)
            results = Eval(dataset, self.input_settings, evaluation.runners)
            evalDF = pd.DataFrame(results)
            #print(evalDF.head())

            outDir = str(self.input_settings.datadir.joinpath(dataset['name'])).replace("inputs/","outputs/")
            out_file = "%s/eval.csv" % (outDir)
            header = True
            if os.path.isfile(out_file):
                print("appending to %s" % (out_file))
                header = False
                # make sure we don't duplicate any rows
                df = pd.read_csv(out_file, header = 0, index_col = 'params')
                evalDF = evalDF[~evalDF.isin(df)].dropna()
            else:
                print("writing to %s" % (out_file))

            with open(out_file, 'a') as out:
                # lock it to avoid scripts trying to write at the same time
                fcntl.flock(out, fcntl.LOCK_EX)
                evalDF.to_csv(out, header=header, index_label='params')
                fcntl.flock(out, fcntl.LOCK_UN)


class ConfigParser(object):
    '''
    Define static methods for parsing a config file that sets a large number
    of parameters for the pipeline
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
            print(algorithm)
            combos = [dict(zip(algorithm['params'], val))
                for val in itertools.product(
                    *(algorithm['params'][param]
                        for param in algorithm['params']))]
            for combo in combos:
                algorithms.append([algorithm['name'],combo])


        return algorithms

    @staticmethod
    def __parse_output_settings(output_settings_map):
        output_dir = Path(output_settings_map['output_dir'])
        return OutputSettings(output_dir)

def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('--config', default='config.yaml',
        help='Configuration file')
    parser.add_argument('--eval-only', action='store_true', default=False,
        help='Skip running the methods and only evaluate')

    return parser


def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts

# In[ ]:


if __name__ == "__main__":
    opts = parse_arguments()
    with open(opts.config, 'r') as conf:
        evaluation = ConfigParser.parse(conf)
    print("%d total runners" % (len(evaluation.runners)))
    #if len(evaluation.runners) > 500 and opts.eval_only is False:
    #    sys.exit("Too many runners. quitting")
    print('Started running methods')

    for idx in range(len(evaluation.runners)):
        evaluation.runners[idx].generateInputs()

    if opts.eval_only is False:
        for idx in range(len(evaluation.runners)):
            status = evaluation.runners[idx].run()
            if status != 'already_exists':
                evaluation.runners[idx].parseOutput()
    print("Finished running")

    evaluation.evaluate_runners()

    print('Evaluation complete')

