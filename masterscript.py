
# Quick script to run/test the algorithms
#print("Importing libraries")

import yaml
#import argparse
from optparse import OptionParser,OptionGroup
from collections import defaultdict
import os
import sys
#from tqdm import tqdm
import itertools
import time
#import networkx as nx
#from scipy import sparse
#import numpy as np
from src.utils import baobab_utils

# import local packages
#import bench


def main(config_map, **kwargs):
    #for alg in config_map['input_settings']['algorithms']:
    #    print(alg)
    #    if alg['name'] == 'GRISLI':
    #        print(alg['params']['L'][1])
    algs = config_map['input_settings']['algorithms']
    print(algs)
    print("algs: %s" % (', '.join(a['name'] for a in algs)))
    #if len(algs) == 1 and algs[0]['name'] == 'SCINGE':
    if len(algs) != 1 or algs[0]['name'] != 'SCINGE':
        print("Currently only setup to run SCINGE. Quitting")
    else:
        # split the SCINGE options into multiple yaml files to combine them for baobab
        params = algs[0]['params']
        params2 = params.copy()
        yaml_base = kwargs['config'].replace('.yaml','')
        os.makedirs(yaml_base, exist_ok=True)
        idx = 1
        for dtnl in params2['dT_num_lags']:
            params['dT_num_lags'] = [dtnl]
            for kw in params2['kernel_width']:
                params['kernel_width'] = [kw]
                for pz in params2['lambda']:
                    params['lambda'] = [pz]

                    yaml_file = "%s/%d.yaml" % (yaml_base, idx)
                    print("writing to %s" % (yaml_file))
                    with open(yaml_file, 'w') as out:
                        yaml.dump(config_map, out, default_flow_style=False)

                    # now make a qsub file and submit
                    qsub_file = os.path.abspath("%s/%d.qsub" % (yaml_base, idx))
                    name = "scinge-%s-%d" % (config_map['input_settings']['datasets'][0]['name'], idx)
                    python = "/data/jeff-law/tools/anaconda3/bin/python"
                    command = "%s -u bench.py --config %s" % (python, os.path.abspath(yaml_file))
                    jobs = ["cd %s" % (os.getcwd()), command]
                    submit = not kwargs['test_run']
                    baobab_utils.writeQsubFile(
                        jobs, qsub_file, name=name, submit=submit,
                        nodes=1, ppn=1, walltime='100:00:00')
                    if kwargs['test_run']:
                        print(qsub_file)
                        sys.exit()

                    # wait for a bit after submitting the first job so it has time to submit the rest
                    if idx == 1:
                        print("waiting for 20 sec after submitting the first job to make sure inputs are setup")
                        time.sleep(20)
                    idx += 1

    #with open(kwargs['config'], 'r') as conf:
    #    evaluation = bench.ConfigParser.parse(conf)

    #for idx in range(len(evaluation.runners)):
    #    if evaluation.runners[idx].name == 'GRISLI':
    #        print(idx, evaluation.runners[idx].params)


def setup_opts():
    ## Parse command line args.
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)

    # general parameters
    group = OptionGroup(parser, 'Main Options')
    group.add_option('','--config', type='string', default="config-files/config.yaml",
                     help="Configuration file")
    group.add_option('','--test-run', action='store_true', default=False,
                     help="Just print out the first command generated")
    parser.add_option_group(group)

    #group = OptionGroup(parser, 'SCINGE Options')
    #group.add_option('-N','--net-file', type='string',
    #                 help="Network file to use. Default is the version's default network")
    #parser.add_option_group(group)

    return parser


def parse_args(args):
    parser = setup_opts()

    (opts, args) = parser.parse_args(args)
    kwargs = vars(opts)
    #kwargs = validate_opts(opts) 
    with open(opts.config, 'r') as conf:
        config_map = yaml.load(conf)
    #ConfigParser.__parse_input_settings(
    #            config_map['input_settings'])
    #ConfigParser.__parse_output_settings(
    #            config_map['output_settings'])
    #print(config_map)

    return config_map, kwargs


if __name__ == "__main__":
    config_map, kwargs = parse_args(sys.argv)
    #print(config_map)
    #sys.exit()

    main(config_map, **kwargs)
