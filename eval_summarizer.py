#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse

import src.plotCurves as pc
import pandas as pd
import eval as ev


class EvalSummarizer(object):
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



    def find_curves(self):
        '''
        Plot PR and ROC curves for each dataset
        for all the algorithms
        '''
        AUPRCDict = {}
        AUROCDict = {}
        uAUPRCDict = {}
        uAUROCDict = {}
        
        
        for dataset in self.input_settings.datasets:
            AUPRC, AUROC, uAUPRC, uAUROC = pc.EvalCurves(dataset, self.input_settings)
            
            AUPRCDict[dataset['name']] = AUPRC
            AUROCDict[dataset['name']] = AUROC
            uAUPRCDict[dataset['name']] = uAUPRC
            uAUROCDict[dataset['name']] = uAUROC
            
        return AUPRCDict, AUROCDict, uAUPRCDict, uAUROCDict

    def time_runners(self):
        '''
        :return:
        dict of times for all dataset, for all algos
        '''
        TimeDict = dict()

        for dataset in self.input_settings.datasets:
            timevals  = self.get_time(dataset)
            TimeDict[dataset["name"]] = timevals

        return TimeDict

    def get_time(self, dataset):
        '''
        :param dataset:
        :return: dict of algorithms, along with the corresponding time for the dataset
        '''
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/" + dataset["name"] + "/"
        algo_dict = dict()
        algos = self.input_settings.algorithms
        for algo in algos:
            path = outDir+algo[0]+"/time.txt"
            time = self.parse_time(path)
            if time == -1:
                print("skipping time computation for ", algo[0], "on dataset", dataset["name"])
                continue

            algo_dict[algo[0]] = time

        return algo_dict


    def parse_time(self, path):
        """
       gets the user time given by the time command, which is stored along with the output when the algorithm runs,
       in a file called time.txt
       :return:
       float containing the time this object took to run on the dataset
        """
        try:
            with open(path, "r") as f:
                lines = f.readlines()
                line = lines[1]
                time_val = float(line.split()[-1])

        except FileNotFoundError:
            print("time output (time.txt) file not found, setting time value to -1")
            time_val = -1
        except ValueError:
            print("Algorithm running failed, setting time value to -1")
            time_val = -1

        return time_val


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

    evaluation = None

    with open(config_file, 'r') as conf:
        evaluation = ev.ConfigParser.parse(conf)
    print(evaluation)
    print('Post-Run Evaluation Summary started')

    eval_summ = EvalSummarizer(evaluation.input_settings,evaluation.output_settings)
    outDir = str(eval_summ.output_settings.base_dir) + \
            str(eval_summ.input_settings.datadir).split("inputs")[1] + "/"+\
            str(eval_summ.output_settings.output_prefix) + "-"

    AUPRCDict, AUROCDict, uAUPRCDict, uAUROCDict = eval_summ.find_curves()

    TimeDict = eval_summ.time_runners()
    
    pd.DataFrame(AUPRCDict).to_csv(outDir+'AUPRCscores.csv')
    pd.DataFrame(AUROCDict).to_csv(outDir+'AUROCscores.csv')
    pd.DataFrame(uAUPRCDict).to_csv(outDir+'uAUPRCscores.csv')
    pd.DataFrame(uAUROCDict).to_csv(outDir+'uAUROCscores.csv')
    pd.DataFrame(TimeDict).to_csv(outDir+'Timescores.csv')



    print('Evaluation complete')


if __name__ == '__main__':
  main()
