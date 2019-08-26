#!/usr/bin/env python

# Quick script to select the most variable genes
# and subset the ExpressionData.csv and refNetwork.csv

import os, sys
import yaml
import argparse
import pandas as pd
#import run_eval_algs


def main(config_map, **kwargs):
    input_settings = config_map['input_settings']
    datasets = input_settings['datasets']
    input_dir = "%s/%s" % (input_settings['input_dir'], input_settings['dataset_dir'])
    algs = input_settings['algorithms']
    #print(algs)
    # first setup the algorithms to run
    #print(kwargs['alg'])
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
            new_config_file = subset_most_var_genes(input_dir, dataset, yaml_base, **kwargs) 
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


def subset_most_var_genes(input_dir, dataset, yaml_base, **kwargs):
    """
    Returns the path to a new yaml file with the subsetting expression matrix
    """
    # first load the GeneOrdering.csv and ExpressionData.csv 
    name = dataset['name']
    dataset_dir = "%s/%s" % (input_dir, name)
    print("\nWorking on %s" % (dataset_dir))
    expr_file = "%s/%s" % (dataset_dir, dataset['exprData'])
    gene_ordering_file = "%s/%s" % (dataset_dir, kwargs['gene_order_file'])

    #postfix = kwargs['pval_cutoff'] if kwargs['pval_cutoff'] is not None else kwargs['num_genes']
    postfix = ""
    postfix += '_'+str(kwargs['pval_cutoff']) if kwargs['pval_cutoff'] is not None else ""
    postfix += "_BF" if kwargs['bf_corr'] else ""
    postfix += '_'+str(kwargs['num_genes']) if kwargs['num_genes'] is not None else ""
    postfix += "_varsort" if kwargs['sort_by_variance'] else ""
    postfix += "_tfs" if kwargs['include_tfs'] else ""
    # remove the leading underscore
    postfix = postfix[1:]
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
    expr_df = pd.read_csv(expr_file, header=0, index_col=0)
    print("\treading %s" % (gene_ordering_file))
    gene_df = pd.read_csv(gene_ordering_file, header=0, index_col=0)
    #print(expr_df.head())
    #print(gene_df.head())
    if kwargs['include_tfs'] is not None:
        print("\treading %s" % (kwargs['include_tfs']))
        tfs_df = pd.read_csv(kwargs['include_tfs'], header=0)
        tfs = tfs_df[tfs_df.columns[0]]
        num_total_tfs = len(tfs)
        # limit the tfs to those present in the expression file
        tfs = tfs[tfs.isin(expr_df.index)]
        print("\t\t%d/%d TFs present in ExpressionData" % (
            len(tfs), num_total_tfs))

    # make sure the variable genes are in the Expression Data
    expr_variable_genes = set(gene_df.index.values) & set(expr_df.index.values)
    extra_genes = set(gene_df.index.values) - set(expr_df.index.values)
    if len(extra_genes) != 0:
        print("WARNING: %d variable genes are not in ExpressionData.csv:" % (len(extra_genes)))
        print(extra_genes)
        gene_df = gene_df.loc[expr_variable_genes]
    # limit the Expression matrix to the variable genes
    pval_col = gene_df.columns[0]
    # make sure its sorted by the pvalue column
    gene_df.sort_values(by=pval_col, inplace=True)
    variable_genes = gene_df.index.values

    # now figure out the genes to subset
    if kwargs['pval_cutoff']:
        pval_cutoff = kwargs['pval_cutoff']
        if kwargs['bf_corr']:
            # divide the pvalue by the # of genes to get the BF-corrected pvalue cutoff
            pval_cutoff = kwargs['pval_cutoff'] / float(len(gene_df.index))
            print("Using the BF-corrected p-value cutoff of %s (%s / %s genes)" % (
                pval_cutoff, kwargs['pval_cutoff'], len(gene_df.index)))

        variable_genes = gene_df[gene_df[pval_col] < pval_cutoff].index.values
        print("\t%d genes pass pval_cutoff of %s" % (len(variable_genes), pval_cutoff))
        print("Before using pValue cut-off num rows: ", gene_df.shape[0] )
        gene_df = gene_df.filter(items = variable_genes, axis='index')
        print("After using pValue cut-off num rows: ", gene_df.shape[0] )
        
    variable_genes = []
    if kwargs['include_tfs'] is not None:
        # include the TFs that pass the p-val cutoff
        tfs = tfs[tfs.isin(gene_df.index)] 
        if kwargs['pval_cutoff']:
            print("\tincluding %d TFs that pass the pval cutoff" % (len(tfs)))
        else:
            print("\tincluding %d TFs" % (len(tfs)))
        variable_tfs = set(tfs)
        gene_df.drop(labels = variable_tfs, axis='index', inplace = True)
    else:
        variable_tfs = set()
    if kwargs['num_genes'] is not None:
        if kwargs['num_genes'] > len(gene_df):
            variable_genes_new = gene_df.index
            pass
        else:
            # they should already be ordered by p-value, so just take them
            if kwargs['sort_by_variance']:
                print("\tsorting by variance")
                if len(gene_df.columns) < 2:
                    print("ERROR: no variance column found. Should be 3rd column. Quitting")
                    sys.exit()
                var_col = gene_df.columns[1]
                print("\tusing the column '%s' as the variance columns" % (var_col))
                # the third column is the variance. Sort by that

                gene_df.sort_values(by=var_col, inplace=True, ascending = False)
             
            variable_genes_new = gene_df.iloc[:kwargs['num_genes']].index.values
        variable_genes = set(variable_genes_new) | set(variable_tfs)

    print("\trestricting to %d genes" % (len(variable_genes)))
    expr_df = expr_df.loc[variable_genes]
    print("\tNew shape of Expression Data %d x %d" % (expr_df.shape[0],expr_df.shape[1]))

    print("\twriting %s" % (new_expr_file))
    expr_df.to_csv(new_expr_file, header=True)
    new_pt_file = "%s/%s" % (new_dataset_dir, dataset['cellData'])
    if not os.path.isfile(new_pt_file):
        # also add a relative symlink to the pseudotime file so its not copied each time
        cwd = os.getcwd()
        os.chdir(new_dataset_dir)
        os.symlink("../"+dataset['cellData'], dataset['cellData'])
        os.chdir(cwd)

    ##### I decided to do this in a separate script since that is 
    ##### really a post-processing step and we want to evaluate multiple networks
    # Also subset the refNetwork.csv file if its there
    #trueEdges_file = "%s/%s" % (dataset_dir, dataset['trueEdges'])
    #if os.path.isfile(trueEdges_file):
    #    print("reading %s" % (trueEdges_file))
    #    trueEdges_df = pd.read_csv(trueEdges_file, header= 0, index_col=None)
    #    print(trueEdges_df.head())

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
    parser.add_argument('--bf-corr', action="store_true", default=False,
        help='Option to correct the p-value using BF correction')
    # TODO specify multiple?
    parser.add_argument('--num-genes', type=int, 
        help='Number of genes to subset. Default: 100')
    parser.add_argument('--sort-by-variance', action="store_true", default=False,
        help='After the p-value cutoff, take the top --num-genes, sorting by the variance. Should be 3rd column in the gene order file')
    parser.add_argument('--include-tfs', type=str,
        help='CSV containing the TFs in the first column. Will include all of them, regardless of the other options set here')
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
