
# Quick script to run/test the algorithms
#print("Importing libraries")

import yaml
#import argparse
from optparse import OptionParser,OptionGroup
from collections import defaultdict
import os
import sys
import subprocess
#from tqdm import tqdm
import itertools
import time
#import networkx as nx
#from scipy import sparse
#import numpy as np
from src.utils import baobab_utils

# import local packages
#import bench
alg_params = {
    "SCODE": "D",
    "SINCERITIES": "nBins",
    "LEAP": "maxLag",
    "SCRIBE": "delay",
    }


def main(config_map, **kwargs):
    #for alg in config_map['input_settings']['algorithms']:
    #    print(alg)
    #    if alg['name'] == 'GRISLI':
    #        print(alg['params']['L'][1])
    algs = config_map['input_settings']['algorithms']
    config_map = config_map.copy()
    print(algs)
    #print("algs: %s" % (', '.join(a['name'] for a in algs)))
    # if there aren't any algs specified by the command line (i.e., kwargs),
    # then use whatever is in the config file
    if kwargs['alg'] is None:
        kwargs['alg'] = [a['name'].lower() for a in algs]
        print("\nNo algs were specified. Using the algorithms in the yaml file:")
        print(str(kwargs['alg']) + '\n')
    else:
        # make the alg names lower so capitalization won't make a difference
        kwargs['alg'] = [a.lower() for a in kwargs['alg']]

    # setup the config file, and then 
    for i, alg in enumerate(algs):
        if alg['name'].lower() in kwargs['alg']:
            print('Running %s' % (alg))
        else:
            continue
        # make a directory for this config, and then put a config file for each method inside
        yaml_base = kwargs['config'].replace('.yaml','')
        os.makedirs(yaml_base, exist_ok=True)
        alg['params']['should_run'] = [True]
        # only write the current alg in this yaml file
        config_map['input_settings']['algorithms'] = [alg]

        if alg['name'] in ['SCINGE', 'GRISLI'] and kwargs['qsub'] is True:
            if alg['name'] == 'SCINGE':
                run_scinge(alg, config_map, yaml_base, **kwargs)
            elif alg['name'] == 'GRISLI':
                run_grisli(alg, config_map, yaml_base, **kwargs)
        else:
            yaml_file = "%s/%s.yaml" % (yaml_base, alg['name'])
            cmd_file = os.path.abspath("%s/%s.sh" % (yaml_base, alg['name']))
            log_file = os.path.abspath("%s/%s.log" % (yaml_base, alg['name']))
            params = alg['params']
            if kwargs['param_per_job']:
                params_list = []
                param_name = alg_params[alg['name']]
                # start a screen session for each parameter
                for p in params[param_name]:
                    params2 = params.copy()
                    params2[param_name] = [p]
                    params_list.append(params2)
            else:
                params_list = [params]
            for curr_params in params_list:
                print(curr_params)
                alg['params'] = curr_params
                write_yaml_file(yaml_file, config_map)
                # now run it. Submit it to screen
                name = "%s-%s" % (alg['name'], config_map['input_settings']['datasets'][0]['name'].split('/')[-1])
                command = "python -u bench.py --config %s %s >> %s 2>&1" % (
                    os.path.abspath(yaml_file), 
                    "--eval-only" if kwargs.get("eval_only") is True else "",
                    log_file,)
                jobs = ["cd %s" % (os.getcwd()), command]
                # write the bash file
                write_bash_file(cmd_file, jobs)
                submit_to_screen(cmd_file, name, log_file, **kwargs)
                # wait a few seconds to make sure the correct yaml file is loaded
                if not kwargs['eval_only']:
                    time.sleep(3)


def write_bash_file(cmd_file, jobs):
    print("\twriting to %s" % (cmd_file))
    with open(cmd_file, 'w') as out:
        out.write('echo "Job Started at: `date`"\n')
        # write each job, as well as an echo (print) statement of the job to be run to the qsub file
        out.write('\n'.join(['echo """%s"""\n%s' % (cmd, cmd) for cmd in jobs]) + '\n')
        out.write('echo "Job Ended at: `date`"\n')


def submit_to_screen(cmd_file, name, log_file, **kwargs):
    print("\tsubmitting '%s' %s to screen" % (name, cmd_file))
    cmd = "screen -S %s -d -m /bin/sh -c \"bash %s >> %s 2>&1\"" % (name, cmd_file, log_file)
    print(cmd+'\n')
    if kwargs['test_run']:
        return
    else:
        subprocess.check_call(cmd, shell=True)


def write_yaml_file(yaml_file, config_map):
    print("\twriting to %s" % (yaml_file))
    with open(yaml_file, 'w') as out:
        yaml.dump(config_map, out, default_flow_style=False)


def run_grisli(alg_settings, config_map, yaml_base, **kwargs):
    alg = alg_settings
    params = alg['params']
    params2 = params.copy()
    idx = 1
    out_dir = "%s/grisli/" % (yaml_base)
    os.makedirs(out_dir, exist_ok=True)
    for a in params2['alphaMin']:
        # this params variable is a pointer to the params inside of config_map
        params['alphaMin'] = [a]

        yaml_file = "%s/%d.yaml" % (out_dir, idx)
        write_yaml_file(yaml_file, config_map)
        # now make a qsub file and submit
        qsub_file = os.path.abspath("%s/%d.qsub" % (out_dir, idx))
        name = "grisli-%s-%d" % (config_map['input_settings']['datasets'][0]['name'], idx)
        python = "/data/jeff-law/tools/anaconda3/bin/python"
        command = "%s -u bench.py --config %s %s" % (
            python, os.path.abspath(yaml_file),
            "--eval-only" if kwargs.get("eval_only") is True else "")
        jobs = ["cd %s" % (os.getcwd()), command]
        submit = not kwargs['test_run']
        baobab_utils.writeQsubFile(
            jobs, qsub_file, name=name, submit=submit,
            nodes=1, ppn=1, walltime='100:00:00')
        if kwargs['test_run']:
            print(qsub_file)
            sys.exit()

        # wait for a bit after submitting the first job so it has time to submit the rest
        if idx == 1 and kwargs.get("eval_only") is not True:
            print("waiting for 20 sec after submitting the first job to make sure inputs are setup")
            time.sleep(20)
        idx += 1


def run_scinge(alg_settings, config_map, yaml_base, **kwargs):
    #if len(algs) == 1 and algs[0]['name'] == 'SCINGE':
    #if len(algs) != 1 or algs[0]['name'] != 'SCINGE':
    #    print("Currently only setup to run SCINGE. Quitting")
    # split the SCINGE options into multiple yaml files to combine them for baobab
    alg = alg_settings
    params = alg['params']
    params2 = params.copy()
    idx = 1
    out_dir = "%s/scinge/" % (yaml_base)
    os.makedirs(out_dir, exist_ok=True)
    for dtnl in params2['dT_num_lags']:
        params['dT_num_lags'] = [dtnl]
        for kw in params2['kernel_width']:
            params['kernel_width'] = [kw]
            for pz in params2['lambda']:
                params['lambda'] = [pz]

                yaml_file = "%s/%d.yaml" % (out_dir, idx)
                write_yaml_file(yaml_file, config_map)
                # now make a qsub file and submit
                qsub_file = os.path.abspath("%s/%d.qsub" % (out_dir, idx))
                name = "scinge-%s-%d" % (config_map['input_settings']['datasets'][0]['name'], idx)
                python = "/data/jeff-law/tools/anaconda3/bin/python"
                command = "%s -u bench.py --config %s %s" % (
                    python, os.path.abspath(yaml_file),
                    "--eval-only" if kwargs.get("eval_only") is True else "")
                jobs = ["cd %s" % (os.getcwd()), command]
                submit = not kwargs['test_run']
                baobab_utils.writeQsubFile(
                    jobs, qsub_file, name=name, submit=submit,
                    nodes=1, ppn=1, walltime='100:00:00')
                if kwargs['test_run']:
                    print(qsub_file)
                    sys.exit()

                # wait for a bit after submitting the first job so it has time to submit the rest
                if idx == 1 and kwargs.get("eval_only") is not True:
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
    group.add_option('','--alg', type='string', action="append", 
                     help="Name of algorithm to run. May specify multiple. Default is whatever is set to true in the config file")
    group.add_option('','--qsub', action='store_true', default=False,
                     help="submit the jobs to a PBS queue with qsub. Currently only setup for GRISLI and SCINGE")
    group.add_option('','--test-run', action='store_true', default=False,
                     help="Just print out the first command generated")
    group.add_option('','--eval-only', action='store_true', default=False,
                     help="Send the --eval-only option to bench.py to only re-evaluate the outputs")
    group.add_option('','--param-per-job', action='store_true', default=False,
                     help="Start a screen session for each parameter of an algorithm")
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
