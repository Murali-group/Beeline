import argparse
import pandas as pd
import BLEval as ev 

def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('config', nargs='+',
        help="List of Configuration file containing list of datasets "
              "algorithms and output specifications.\n")
 
    parser.add_argument('-a', '--auc', action="store_true", default=False,
        help="Compute median of areas under Precision-Recall and ROC curves.\n")
    
    parser.add_argument('-j', '--jaccard', action="store_true", default=False,
      help="Compute median Jaccard index of predicted top-k networks "
      "for each algorithm for a given set of datasets generated "
      "from the same ground truth network.\n")

    parser.add_argument('-r', '--spearman', action="store_true", default=False,
      help="Compute median Spearman Corr. of predicted edges "
      "for each algorithm  for a given set of datasets generated "
      " from the same ground truth network.\n")

    parser.add_argument('-t', '--time', action="store_true", default=False,
      help="Analyze time taken by each algorithm for a.\n")
    
    parser.add_argument('-e', '--epr', action="store_true", default=False,
      help="Compute median early precision.")
    
    parser.add_argument('-s','--sepr', action="store_true", default=False,
      help="Analyze median (signed) early precision for activation and inhibitory edges.")

      
    parser.add_argument(
        '-o', '--output_dir',
        required=True,
        type=str,
        help='Path to Output Directory',
    )
    return parser

def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()
    
    return opts

def output_prefix(config_file):
    with open(config_file, 'r') as conf:
        evalConfig = ev.ConfigParser.parse(conf)
        evalSummarizer = ev.BLEval(
            evalConfig.input_settings,
            evalConfig.output_settings
        )
        return str(evalSummarizer.output_settings.base_dir) + \
            str(evalSummarizer.input_settings.datadir).split("inputs")[1] + "/"+\
            str(evalSummarizer.output_settings.output_prefix) + "-"


def collate_tables(data_dirs, df_suffix, output_dir):
    csv_files = [outDir+df_suffix for outDir in data_dirs]
    csv_dfs = [pd.read_csv(fx, index_col=0) for fx in csv_files]
    merged_df = pd.concat(csv_dfs, axis=1).T
    merged_df.to_csv(output_dir+df_suffix)
    return merged_df


def main():
    opts = parse_arguments()
    print(opts)
    evalOutDirs = [output_prefix(cfg_file) for cfg_file in opts.config] 
    print('\nPost-run evaluation started...')
    output_dir = opts.output_dir + '/'
    df_dict = {}

    # Compute and plot ROC, PRC and report median AUROC, AUPRC    
    if (opts.auc):
        print('\n\nComputing areas under ROC and PR curves...')
        df_dict['AUPR'] = collate_tables(evalOutDirs, 'AUPRC.csv', output_dir)
        df_dict['AUROC'] = collate_tables(evalOutDirs, 'AUROC.csv', output_dir)

    # Compute median time taken
    if (opts.time):
        print('\n\nComputing time taken...')
        df_dict['Time'] = collate_tables(evalOutDirs, 'Times.csv', output_dir)
    
    # Compute early precision
    if (opts.epr):
        print('\n\nComputing early precision values...')
        df_dict['Early Percision'] = collate_tables(evalOutDirs, 'EPr.csv', output_dir)
  
    # Compute Jaccard index    
    if (opts.jaccard):
        print('\n\nComputing Jaccard index...')
        df_dict['Jaccard'] = collate_tables(evalOutDirs, "Jaccard.csv", output_dir)

    # Compute Spearman correlation scores
    if (opts.spearman):
        print('\n\nComputing Spearman\'s correlation...')
        df_dict['Spearman'] = collate_tables(evalOutDirs, "Spearman.csv", output_dir)
 
     # Compute early precision for activation and inhibitory edges
    if (opts.sepr):
        print('\n\nComputing early precision values for activation and inhibitory edges...')
        df_dict['EPR Activation'] = collate_tables(evalOutDirs, "EPr-Activation.csv", output_dir)
        df_dict['EPR Inhibitory'] = collate_tables(evalOutDirs, "EPr-Inhibitory.csv", output_dir)

    output_xlsx = output_dir + 'BEELINE-Collated.xlsx'
    with pd.ExcelWriter(output_xlsx) as writer:
        for dfkey, rdf in df_dict.items():
            rdf.to_excel(writer, sheet_name=dfkey)


if __name__ == "__main__":
    main()
