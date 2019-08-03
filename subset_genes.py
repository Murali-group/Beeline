#!/usr/bin/env python

# Quick script to select the most variable genes
# and subset the ExpressionData.csv and refNetwork.csv

import os
import yaml
import argparse
import pandas as pd
import run_eval_algs


def main(config_map, **kwargs):
    input_settings = config_map['input_settings']
    datasets = input_settings['datasets']
    input_dir = "%s/%s" % (input_settings['input_dir'], input_settings['dataset_dir'])
    algs = input_settings['algorithms']
    #print(algs)
    # first setup the algorithms to run
    print(kwargs['alg'])
    if kwargs['alg'] is not None:
        # make the alg names lower so capitalization won't make a difference
        kwargs['alg'] = [a.lower() for a in kwargs['alg']]
        new_alg_settings = []
        for alg in algs:
            if alg['name'].lower() in kwargs['alg']:
                print('Keeping %s in the new config files' % (alg))
            else:
                continue
            # set 'should_run' to True for the algs specified
            alg['params']['should_run'] = [True]
            new_alg_settings.append(alg)
        input_settings['algorithms'] = new_alg_settings
    config_map = config_map.copy()

    # make a directory for this config, and then put a config file for each method inside
    yaml_base = kwargs['config'].replace('.yaml','')
    os.makedirs(yaml_base, exist_ok=True)
    new_config_files = []

    if kwargs['most_variable_genes'] is True:
        for dataset in datasets:
            new_config_file = subset_most_var_genes(input_dir, dataset, yaml_base) 
            # only write the current dataset in this yaml file
            config_map['input_settings']['datasets'] = [dataset]
            write_yaml_file(new_config_file, config_map) 
            new_config_files.append(new_config_file) 

    # TODO allow running without subsetting genes
    #if len(new_config_files) == 0:
    #    new_config_files.append(kwargs['config'])
    #    new_config_file = "%s/%s.yaml" % (yaml_base, name, postfix)
    #    write_yaml_file(new_config_file, config_map) 

    if kwargs.get('run_algs') is True:
        # after the datasets are setup, run the methods on the created config files
        for config_file in new_config_files:
            print("\nRunning run_eval_algs.py on %s" % (config_file))
            with open(config_file, 'r') as conf:
                config_map = yaml.load(conf)

            # TODO start a screen session to run this 
            run_eval_algs.main(config_map)


def subset_most_var_genes(input_dir, dataset, yaml_base):
    """
    Returns the path to a new yaml file with the subsetting expression matrix
    """
    # first load the GeneOrdering.csv and ExpressionData.csv 
    name = dataset['name']
    dataset_dir = "%s/%s" % (input_dir, name)
    print("\nWorking on %s" % (dataset_dir))
    expr_file = "%s/%s" % (dataset_dir, dataset['exprData'])
    gene_ordering_file = "%s/%s" % (dataset_dir, kwargs['gene_order_file'])

    postfix = kwargs['pval_cutoff'] if kwargs['pval_cutoff'] is not None else kwargs['num_genes']
    #new_expr_file = "%s/%s" % (dataset_dir, dataset['exprData'].replace('.csv', '-%s.csv' % (postfix)))
    # put it in a new directory so the outputs don't overwrite each other
    dataset['name'] += "/genecutoff-%s" % (postfix)
    new_dataset_dir = "%s/%s" % (input_dir, dataset['name'])
    os.makedirs(new_dataset_dir, exist_ok=True)
    new_expr_file = "%s/%s" % (new_dataset_dir, dataset['exprData'].replace('.csv', '-%s.csv' % (postfix)))
    new_config_file = "%s/%s-%s.yaml" % (yaml_base, name, postfix)
    # change the expression data in the config file
    dataset['exprData'] = new_expr_file.split('/')[-1]
    if os.path.isfile(new_expr_file) and kwargs.get('forced') is not True:
        print("\t%s already exists. Use --forced to overwrite" % (new_expr_file))
        return new_config_file

    print("\treading %s" % (expr_file))
    expr_df = pd.read_csv(expr_file, header= 0, index_col=0)
    print("\treading %s" % (gene_ordering_file))
    gene_df = pd.read_csv(gene_ordering_file, header= 0, index_col=0)
    #print(expr_df.head())
    #print(gene_df.head())

    # make sure the variable genes are in the Expression Data
    expr_variable_genes = set(gene_df.index.values) & set(expr_df.index.values)
    extra_genes = set(gene_df.index.values) - set(expr_df.index.values)
    if len(extra_genes) != 0:
        print("WARNING: %d variable genes are not in ExpressionData.csv:" % (len(extra_genes)))
        print(extra_genes)
        gene_df = gene_df.loc[expr_variable_genes]
    # limit the Expression matrix to the variable genes
    if kwargs['pval_cutoff']:
        pval_col = gene_df.columns[0]
        variable_genes = gene_df[gene_df[pval_col] < kwargs['pval_cutoff']].index.values
    elif kwargs['num_genes']:
        variable_genes = gene_df.iloc[:kwargs['num_genes']].index.values

    print("\trestricting to %d genes" % (len(variable_genes)))
    expr_df = expr_df.loc[variable_genes]
    print("\twriting %s" % (new_expr_file))
    expr_df.to_csv(new_expr_file, header=True)
    # also add a relative symlink to the pseudotime file so its not copied each time
    cwd = os.getcwd()
    os.chdir(new_dataset_dir)
    os.symlink("../"+dataset['cellData'], dataset['cellData'])
    os.chdir(cwd)

    # also subset the refNetwork.csv file if its there
    trueEdges_file = "%s/%s" % (dataset_dir, dataset['trueEdges'])
    if os.path.isfile(trueEdges_file):
        print("reading %s" % (trueEdges_file))
        trueEdges_df = pd.read_csv(trueEdges_file, header= 0, index_col=None)
        print(trueEdges_df.head())

    # now write the config file
    return new_config_file


def write_yaml_file(yaml_file, config_map):
    print("\twriting to %s" % (yaml_file))
    with open(yaml_file, 'w') as out:
        yaml.dump(config_map, out, default_flow_style=False)


def setup_parser():
    parser = argparse.ArgumentParser(
        description='Script for setting up various experiments ')

    parser.add_argument('--config', default='config.yaml', required=True,
        help='Configuration file')
    parser.add_argument('--run-algs', action="store_true", default=False,
        help='Run the methods using the generated config file')
    parser.add_argument('--alg', action="append", 
        help="Name of algorithm to run. May specify multiple. Default is whatever is set to true in the config file")

    # most variable genes options
    parser.add_argument('--most-variable-genes', '-V', action="store_true", default=False,
        help='Select the most variable genes and subset the Expression Data.csv and refNetwork.csv to those genes')
    parser.add_argument('--gene-order-file', type=str, default="GeneOrdering.csv",
        help='Name of CSV file with the ascending ordering value in the second column. ' +
        'Should be the same for each dataset. Suggested: GeneOrdering.csv.')
    # TODO specify multiple?
    parser.add_argument('--pval-cutoff', type=float, 
        help='Cutoff of the pvalue to select genes')
    # TODO specify multiple?
    parser.add_argument('--num-genes', type=int, default=100,
        help='Number of genes to subset. Default: 100')
    parser.add_argument('--forced', action="store_true", default=False,
        help='Overwrite the ExpressionData.csv file if it already exists.')

    return parser


if __name__ == "__main__":
    parser = setup_parser()
    opts = parser.parse_args()
    kwargs = vars(opts)
    config_file = opts.config
    with open(config_file, 'r') as conf:
        config_map = yaml.load(conf)

    main(config_map, **kwargs)
