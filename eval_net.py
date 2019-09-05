#!/usr/bin/env python

# Script to limit the ground truth network to the genes which are in the 
# expression data file, and evaluate


import os
import yaml
import argparse
import pandas as pd
#import run_eval_algs
import BLEvalAggregator as BLeval


def main(config_map, opts):
    config_map = config_map.copy()
    input_settings = config_map['input_settings']
    out_settings = config_map['output_settings']
    datasets = input_settings['datasets']
    input_dir = "%s/%s" % (input_settings['input_dir'], input_settings['dataset_dir'])
    algs = input_settings['algorithms']
    if opts.alg is not None:
        # make the alg names lower so capitalization won't make a difference
        opts.algs = [a.lower() for a in opts.alg]
        new_alg_settings = []
        for alg in opts.alg:
            # set 'should_run' to True for the algs specified
            algdict = {'name': alg, 'params': {'should_run': [True]}}
            new_alg_settings.append(algdict)
        input_settings['algorithms'] = new_alg_settings

    for dataset in datasets:
        # first load ExpressionData.csv 
        name = dataset['name']
        dataset_dir = "%s/%s" % (input_dir, name)
        print("\nWorking on %s" % (dataset_dir))
        expr_file = "%s/%s" % (dataset_dir, dataset['exprData'])
        print("\treading %s" % (expr_file))
        expr_df = pd.read_csv(expr_file, header= 0, index_col=0)

        # now load the network file
        net_file = "%s/%s" % (dataset_dir, dataset['trueEdges'])
        print("\treading %s" % (opts.ref_net_file))
        net_df = pd.read_csv(opts.ref_net_file, header=0)
        net_df.columns = ["Gene1","Gene2"] + list(net_df.columns[2:])
        net_tfs = net_df['Gene1'].values
        num_tfs, num_targets = net_df[['Gene1','Gene2']].nunique()
        print("\t%d TFs, %d targets, %d edges"  % (num_tfs, num_targets, len(net_df)))

        expr_genes = set(expr_df.index.values)

        net_df = net_df[(net_df['Gene1'].isin(expr_genes) & net_df['Gene2'].isin(expr_genes))]
        
                
        if (opts.tfs_only):
            net_df = net_df[(net_df['Gene2'].isin(net_df['Gene1'].unique()))]            
            
            
        if len(net_df) == 0:
            print("No matching node names found. Please make sure the same namespace is used.")
            print("\tExample expr node: %s" % (list(expr_genes)[0]))
            print("\tExample net node: %s" % (net_tfs[0]))
        else:
            print("After limitting to the %d genes with expression values:" % (len(expr_genes)))
            num_tfs, num_targets = net_df[['Gene1','Gene2']].nunique()
            print("\t# TFs\t# targets\t# edges")
            print("\t%s\t%s\t%d"  % (num_tfs, num_targets, len(net_df)))

            outFile = open('stats-new.txt','a')
            
            if opts.stats_only:
                
                outFile.write("\t%s\t%s\t%d\t%s\t%s\t%s\n"  % (num_tfs, num_targets, len(net_df), dataset['name'].split('/genecutoff-0.01_BF_')[0],dataset['name'].split('/genecutoff-0.01_BF_')[1], opts.ref_net_file))
                continue
            # and write it to a file

            print("\nwriting %s" % (net_file))
            net_df.to_csv(net_file, index=False)

            
            # don't need to write the yaml file
            # add an option to write it?
            # can simply pass it to BLEvalAggregator.py
            print("Running BLEvalAggregator.py")
            print(opts)
            BLeval.main(config_map, opts)

            # after its done, need to move the evaluation file
            # otherwise it will be overwritten by the next run
            # alternatively we could change the output directory in the config map

            net_name = opts.net_name if opts.net_name is not None else opts.ref_net_file.split('/')[-1].replace('.csv','')
            out_file = input_dir.replace("inputs/","outputs/")+"/outStats/"+str(opts.net_name)+".csv"
            print(out_file)
            all_df = pd.DataFrame()
            #for measure in ["AUPRC", "AUROC", "EPr", "Jaccard", "Times"]:
            for measure in ["AUPRC", "EPr", "Times","trueEC","possibleEC","nodeCnt","TFCnt"]:
                measure_file = "%s/%s-%s.csv" % (
                    input_dir.replace("inputs/","outputs/"), out_settings['output_prefix'], measure)
                
                df = pd.read_csv(measure_file, header=0)
                df.columns = ['algorithm', 'value']
                df['measure'] = measure
                df['dataset'] = dataset['name']
                df['ref_net'] = net_name
                all_df = pd.concat([all_df, df])
                # delete this file
                os.remove(measure_file)

            # now append this to a file
            header = True
            append = True
            if os.path.isfile(out_file):
                #if forced:
                #    append = False
                #    print("writing to %s" % (out_file))
                #else:
                print("appending to %s" % (out_file))
                #header = False
                # make sure we don't duplicate any rows
                df = pd.read_csv(out_file, header = 0)
                all_df = pd.concat([df, all_df])
                # if the new values are already in the df, don't repeat them again
                all_df.drop_duplicates(inplace=True)
                # if the new values are different, overwrite what was in the file with the new results
                all_df.drop_duplicates(subset=["algorithm", "measure", "dataset", "ref_net"], keep='last', inplace=True)
            else:
                print("writing to %s" % (out_file))

            #with open(out_file, 'a' if append else 'w') as out:
            with open(out_file, 'w') as out:
                # lock it to avoid scripts trying to write at the same time
                #fcntl.flock(out, fcntl.LOCK_EX)
                all_df.to_csv(out, header=header, index=False)
                #fcntl.flock(out, fcntl.LOCK_UN)
    #outFile.close()
    print("Finished")


def write_yaml_file(yaml_file, config_map):
    print("\twriting to %s" % (yaml_file))
    with open(yaml_file, 'w') as out:
        yaml.dump(config_map, out, default_flow_style=False)


def setup_parser():
    #parser = argparse.ArgumentParser(
    #    description='Script for setting up various experiments ')
    # also add the BLEval options
    parser = BLeval.get_parser()

    #parser.add_argument('--config', default='config.yaml', required=True,
    #    help='Configuration file')
    #parser.add_argument('--run-algs', action="store_true", default=False,
    #    help='Run the methods using the generated config file')
    parser.add_argument('--alg', action="append", 
        help="Name of algorithm to run. Must match the output file path. May specify multiple. Default is whatever is set to true in the config file")
    
    parser.add_argument('--ref-net-file', type=str, default="GeneOrdering.csv",
        help='Path to the ground truth refNetwork.csv file. A new file will be subset to the genes in the ExpressionData.csv and written.')
    
    parser.add_argument('--net-name', type = str, 
        help='The name to give this network for evaluating. Default is the file name.')
    
    parser.add_argument('--stats-only', action="store_true", default=False,
        help='Only print out the stats of the # edges and such')

    
    parser.add_argument('--tfs-only', action="store_true", default=False,
        help='Only use TFs which have outgoing edges as nodes.')
    ## most variable genes options
    #parser.add_argument('--most-variable-genes', '-V', action="store_true", default=False,
    #    help='Select the most variable genes and subset the Expression Data.csv and refNetwork.csv to those genes')
    #parser.add_argument('--gene-order-file', type=str, default="GeneOrdering.csv",
    #    help='Name of CSV file with the ascending ordering value in the second column. ' +
    #    'Should be the same for each dataset. Suggested: GeneOrdering.csv.')
    # TODO specify multiple?
    #parser.add_argument('--pval-cutoff', type=float, 
    #    help='Cutoff of the pvalue to select genes')
    # TODO specify multiple?
    #parser.add_argument('--num-genes', type=int, default=100,
    #    help='Number of genes to subset. Default: 100')
    #parser.add_argument('--forced', action="store_true", default=False,
    #    help='Overwrite the ExpressionData.csv file if it already exists.')

    return parser


if __name__ == "__main__":
    parser = setup_parser()
    opts = parser.parse_args()
    
    # BLEval takes the opts, so keep it as opts
    #kwargs = vars(opts)
    config_file = opts.config
    with open(config_file, 'r') as conf:
        config_map = yaml.load(conf)

    main(config_map, opts)
    
