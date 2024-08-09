#!/usr/bin/env python
# coding: utf-8

# Please refer to https://github.com/theislab/sfaira_tutorials/blob/master/tutorials/data_loaders.ipynb

import argparse
import tensorflow as tf

# local imports
import BLData as dt

def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(description='Download scRNA-seq datasets from Sfaira.')

    # Specify configure file
    parser.add_argument('--config', default='config.yaml',
        help="Configuration file containing list of input setting specifications.\n")

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

    with open(config_file, 'r') as conf:
        sfairaLoader = dt.ConfigParser.parse(conf)
    print(sfairaLoader)
    
    dataSummarizer = dt.SfairaData(sfairaLoader.sfaira_settings)

    print("## Dataset downloads started")
    dataSummarizer.sfairaLoader()
    print('##Dataset downloads complete')
    

if __name__ == '__main__':
  main()
