
# script to write a config yaml file with the best performing parameters for each alg and each dataset

import os, sys
import yaml
import argparse
import pandas as pd
import subprocess
#import bench
import shutil


# mapping from the string I used to the actual parameter in the config file
param_mapping = {
# grisli
    "a": 'alphaMin',
    # scinge
    'l': 'lambda',
    'nl': 'num_lags',
    'kw': 'kernel_width',
    'pz': 'prob_zero_removal',
    'pr': 'prob_remove_samples',
    'nr': 'num_replicates',
    }


def main(config_map, **kwargs):
    input_settings = config_map['input_settings']
    algs = input_settings['algorithms']
    dataset = input_settings['datasets'][0]
    config_map = config_map.copy()
    #print(algs)
    ##print("algs: %s" % (', '.join(a['name'] for a in algs)))
    ## if there aren't any algs specified by the command line (i.e., kwargs),
    ## then use whatever is in the config file
    #if kwargs['alg'] is None:
    #    kwargs['alg'] = [a['name'].lower() for a in algs]
    #    print("\nNo algs were specified. Using the algorithms in the yaml file:")
    #    print(str(kwargs['alg']) + '\n')
    #else:
    #    # make the alg names lower so capitalization won't make a difference
    #    kwargs['alg'] = [a.lower() for a in kwargs['alg']]

    print("reading %s" % (kwargs['best_params_file']))
    best_params_df = pd.read_csv(kwargs['best_params_file'], header=0, index_col=0)
    print(best_params_df.head())

    dataset_df = None
    for dataset_name, df_group in best_params_df.groupby("dataset"):
        print(dataset_name)
        if dataset_name in dataset['name']:
            print("Using '%s'" % (dataset_name))
            dataset_df = df_group 
            break
    if dataset_df is None:
        print("ERROR: dataset not found in the best params file %s" % (kwargs['best_params_file']))
        sys.exit("Quitting")
    print(dataset_df.head())

    new_algs = []
    # for each algorithm, get the best parameters
    for alg in algs:
        name = alg['name']
        if name not in dataset_df.index:
            continue
        # extract the parameters from the params string
        params_str = dataset_df.loc[name]['parameter']
        print(name, params_str)
        # this is kind of hacky, but it works for now
        if name == "SCRIBE":
            p = "delay"
            params = {p: params_str.replace(p,'')}
        else:
            params = get_param_list(params_str)
        # map the parameter names to the names if they're there
        params = {param_mapping.get(p,p): val for p, val in params.items()}
        #print(params)
        for param, val in params.items():
            alg['params'][param] = [val]
        if name in ["SCINGE", "SINGE"]:
            alg['params'].pop("dT_num_lags")
        alg['params']['should_run'] = [True]
        # the main branch still uses z instead of D
        if name == "SCODE":
            alg['params']['z'] = alg['params']['D']
            del alg['params']['D']
        print(alg)
        new_algs.append(alg) 

    config_map['input_settings']['algorithms'] = new_algs
    # get the output file for this alg, and copy it to the output dir
    #evaluation = bench.ConfigParser.parse(config_map)
    ##print(config_map['input_settings'])
    #for idx in range(len(evaluation.runners)):
    #    #new_out_dir = "%s/%s/%s/"
    #    runner = evaluation.runners[idx]
    #    runner.setupParams()
        #print(runner.final_ranked_edges)
        #dataset = str(runner.final_ranked_edges).split('/')[2]
        #genes = str(runner.final_ranked_edges).split('/')[3]
        #new_out_file = "%s/%s/%s/%s/rankedEdges.csv" % (kwargs['copy_out_files'], dataset, genes, runner.name)
        #print(new_out_file)
        #os.makedirs(os.path.dirname(new_out_file), exist_ok=True)
        #shutil.copy(runner.final_ranked_edges,new_out_file)

    # make a directory for this config, and then put a config file for each method inside
    yaml_base = kwargs['config'].replace('.yaml','')
    #os.makedirs(yaml_base, exist_ok=True)
    yaml_file = "%s-best-params.yaml" % (yaml_base)
    write_yaml_file(yaml_file, config_map)


def get_param_list(params_str):
    params = {}
#     print(df.head())
    for i, param_val in enumerate(params_str.split('-')):
        # the parameter name could be multiple letters. Find the first index of a number
        num_idx = min([j for j in range(len(param_val)) if param_val[j].isdigit()])
        p = param_val[:num_idx]
        val = param_val[num_idx:]
        params[p] = val
    return params


def write_yaml_file(yaml_file, config_map):
    print("\twriting to %s" % (yaml_file))
    with open(yaml_file, 'w') as out:
        yaml.dump(config_map, out, default_flow_style=False)


def setup_parser():
    parser = argparse.ArgumentParser(
        description='Script for setting up various experiments ')

    parser.add_argument('--config', default='config.yaml', required=True,
        help='Configuration file')
    #parser.add_argument('--alg', action="append", 
    #    help="Name of algorithm to run. May specify multiple. Default is whatever is set to true in the config file")

    # most variable genes options
    parser.add_argument('--best-params-file', type=str, required=True,
        help='path to CSV with the these columns: alg, dataset, params.')
    parser.add_argument('--copy-out-files', 
        help="Copy the output files to another directory")

    return parser


if __name__ == "__main__":
    parser = setup_parser()
    opts = parser.parse_args()
    kwargs = vars(opts)
    config_file = opts.config
    with open(config_file, 'r') as conf:
        config_map = yaml.load(conf)

    main(config_map, **kwargs)

