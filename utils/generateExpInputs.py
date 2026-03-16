from __future__ import print_function

import argparse
import copy
import sys
from pathlib import Path

import pandas as pd
import yaml


def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Generate experimental scRNA-seq inputs for BEELINE.',
        epilog=(
            'Example usage to generate dataset with all TFs and 500 genes (TFs+500): '
            'python generateExpInputs.py --outputDirectory=outputs --dataset-id=hESC-500 '
            '-e=ExpressionData.csv -g=GeneOrdering.csv -P=PseudoTime.csv '
            '-f=STRING-network.csv -i=human-tfs.csv -p=0.01 -n=500'
        )
    )

    parser.add_argument('-e', '--expFile', type=str,
                        default='ExpressionData.csv',
                        help='Path to expression data file. Required.\n')

    parser.add_argument('-g', '--geneOrderingFile', type=str,
                        default='GeneOrdering.csv',
                        help='Path to gene ordering file. Required.\n')

    parser.add_argument('-P', '--pseudoTimeFile', type=str,
                        default='PseudoTime.csv',
                        help='Path to pseudotime file. Copied to output as PseudoTime.csv. Required.\n')

    # default=None so omitting -f skips network processing entirely.
    parser.add_argument('-f', '--netFile', type=str,
                        default=None,
                        help='Path to network file to filter and produce GroundTruthNetwork.csv. Optional.\n')

    parser.add_argument('-i', '--TFFile', type=str,
                        default='human-tfs.csv',
                        help='Path to file containing list of TFs. Required.\n')

    parser.add_argument('-p', '--pVal', type=float, default=0.01,
                        help='p-value cutoff. Default = 0.01')

    parser.add_argument('-n', '--numGenes', type=int, default=500,
                        help='Number of non-TF genes to include. Default=500.\n')

    parser.add_argument('--outputDirectory', type=str, default='outputs',
                        help='Base output directory. Outputs are written to outputDirectory/dataset-id/.\n')

    parser.add_argument('--dataset-id', dest='dataset_id', type=str, default='dataset',
                        help='Identifier for this dataset. Used as the output subdirectory name '
                             'and as the dataset_id column in dataset_stats.csv.\n')

    # Boolean flag pairs: pass the positive flag to enable, --no-X to disable.
    parser.add_argument('-c', '--BFcorr', dest='BFcorr', action='store_true',
                        help='Perform Bonferroni correction on the p-value cutoff (default).\n')
    parser.add_argument('--no-BFcorr', dest='BFcorr', action='store_false',
                        help='Disable Bonferroni correction.\n')
    parser.set_defaults(BFcorr=True)

    parser.add_argument('-t', '--TFs', dest='TFs', action='store_true',
                        help='Add all significantly varying TFs to the gene set (default).\n')
    parser.add_argument('--no-TFs', dest='TFs', action='store_false',
                        help='Do not force-include TFs.\n')
    parser.set_defaults(TFs=True)

    parser.add_argument('-s', '--sort-variance', dest='sortByVariance', action='store_true',
                        help='Select top-n non-TF genes by variance (default).\n')
    parser.add_argument('--no-sort-variance', dest='sortByVariance', action='store_false',
                        help='Select top-n non-TF genes by p-value order instead of variance.\n')
    parser.set_defaults(sortByVariance=True)

    parser.add_argument('-y', '--config', type=str, default=None,
                        help='Path to a YAML config file. When provided, all other arguments '
                             'are sourced from the file; any CLI args supplied alongside -y '
                             'are ignored (a warning is printed).\n')

    return parser


def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()
    return opts


_YAML_DATASET_REQUIRED_KEYS = {
    'dataset_id', 'expFile', 'geneOrderingFile', 'pseudoTimeFile', 'TFFile',
    'pVal', 'BFcorr', 'TFs', 'numGenes', 'sortByVariance',
}
# netFile is optional per dataset (may be absent or null).


def load_yaml_config(path: str, base_opts: argparse.Namespace) -> list:
    '''
    Load a YAML config file and return one argparse.Namespace per dataset entry.

    The YAML must have a top-level 'settings' block with 'outputDirectory' and
    a top-level 'datasets' list. Each dataset entry must contain all keys in
    _YAML_DATASET_REQUIRED_KEYS; 'netFile' is optional (defaults to None).

    The outPrefix for each dataset is derived as outputDirectory/dataset_id.

    Warns if any CLI arguments other than -y/--config were supplied alongside
    the YAML path, because the YAML takes full precedence over them.

    :param path:      str. Path to the YAML config file.
    :param base_opts: argparse.Namespace. Parsed CLI arguments; used only as a
                      template — fields are overwritten per dataset from the YAML.
    :return:          list of argparse.Namespace, one entry per dataset.
    '''
    if not isinstance(path, str):
        raise TypeError(f"path must be str, got {type(path)}")
    if not isinstance(base_opts, argparse.Namespace):
        raise TypeError(f"base_opts must be argparse.Namespace, got {type(base_opts)}")

    # Warn if the user also passed non-config CLI args — YAML wins.
    extra_args = [
        a for a in sys.argv[1:]
        if a not in ('-y', '--config') and not a.endswith('.yaml') and not a.endswith('.yml')
    ]
    if extra_args:
        print(
            "WARNING: CLI arguments %s were provided alongside -y/--config. "
            "The YAML config file takes precedence; these arguments will be ignored." % extra_args
        )

    with open(path, 'r') as f:
        cfg = yaml.safe_load(f)

    if 'settings' not in cfg or 'outputDirectory' not in cfg.get('settings', {}):
        print("ERROR: YAML config must have a 'settings.outputDirectory' field.")
        sys.exit(1)
    if 'datasets' not in cfg or not isinstance(cfg['datasets'], list):
        print("ERROR: YAML config must have a 'datasets' list.")
        sys.exit(1)

    output_dir = cfg['settings']['outputDirectory']
    all_opts = []

    for i, ds in enumerate(cfg['datasets']):
        missing = _YAML_DATASET_REQUIRED_KEYS - set(ds.keys())
        if missing:
            print("ERROR: dataset entry %d is missing required keys: %s" % (i, sorted(missing)))
            sys.exit(1)

        opts = copy.copy(base_opts)
        opts.dataset_id       = ds['dataset_id']
        opts.outputDirectory  = output_dir
        opts.expFile          = ds['expFile']
        opts.geneOrderingFile = ds['geneOrderingFile']
        opts.pseudoTimeFile   = ds['pseudoTimeFile']
        opts.netFile          = ds.get('netFile', None)
        opts.TFFile           = ds['TFFile']
        opts.pVal             = float(ds['pVal'])
        opts.BFcorr           = bool(ds['BFcorr'])
        opts.TFs              = bool(ds['TFs'])
        opts.numGenes         = int(ds['numGenes'])
        opts.sortByVariance   = bool(ds['sortByVariance'])
        opts.outPrefix        = str(Path(output_dir) / ds['dataset_id'])
        all_opts.append(opts)

    return all_opts


def filter_genes_by_pvalue(
    gene_df: pd.DataFrame,
    pval_cutoff: float,
    bf_corr: bool
) -> pd.DataFrame:
    '''
    Filter a gene ordering DataFrame to retain only genes whose p-value falls
    below the given cutoff. Optionally applies Bonferroni correction to the
    cutoff before filtering.

    The p-value column is assumed to be the first data column of gene_df (i.e.,
    gene_df.columns[0]). The DataFrame is expected to already be sorted by this
    column in ascending order before calling this function.

    When bf_corr is True the cutoff is divided by the total number of genes in
    gene_df, so the effective per-gene alpha is pval_cutoff / len(gene_df).
    This controls the family-wise error rate across all tests.

    If pval_cutoff is 0 no filtering is applied and gene_df is returned
    unchanged.

    :param gene_df:     pd.DataFrame indexed by gene name. First column must
                        contain p-values (float). Additional columns are
                        preserved unchanged.
    :param pval_cutoff: float. Nominal p-value threshold. 0 disables filtering.
    :param bf_corr:     bool. When True, divide pval_cutoff by len(gene_df)
                        before filtering (Bonferroni correction).
    :return:            pd.DataFrame with the same columns as gene_df, containing
                        only rows whose p-value is strictly less than the
                        (possibly corrected) cutoff.
    '''
    if pval_cutoff == 0:
        return gene_df

    pval_col = gene_df.columns[0]

    if bf_corr:
        # Divide the nominal alpha by the number of genes to get the
        # per-comparison threshold that controls family-wise error rate.
        corrected_cutoff = pval_cutoff / float(len(gene_df.index))
        print("\nUsing the BF-corrected p-value cutoff of %s (%s / %s genes)" % (
            corrected_cutoff, pval_cutoff, len(gene_df.index)))
        pval_cutoff = corrected_cutoff

    filtered = gene_df[gene_df[pval_col] < pval_cutoff]
    print("\n%d genes pass pval_cutoff of %s" % (len(filtered), pval_cutoff))
    print("\nAfter using pValue cut-off num rows: ", filtered.shape[0])
    return filtered


def select_gene_set(
    gene_df: pd.DataFrame,
    variable_tfs: set,
    num_genes: int,
    sort_by_variance: bool
) -> list:
    '''
    Build the final gene set by combining a pool of non-TF variable genes with
    a pre-determined set of TF genes.

    Non-TF genes are drawn from gene_df after TF rows have been removed. Up to
    num_genes non-TF genes are selected; when sort_by_variance is True the top
    candidates are those with the highest variance (second data column of
    gene_df), otherwise they are taken in p-value order (i.e. the current row
    order of gene_df, assumed already sorted ascending by p-value).

    When num_genes is 0, no non-TF genes are included and the returned list
    contains only TFs.

    When num_genes exceeds the number of available non-TF genes in gene_df, all
    available genes are included without truncation.

    :param gene_df:          pd.DataFrame indexed by gene name. Must contain at
                             least one column (p-value). When sort_by_variance is
                             True a second column (variance, float) is also
                             required; the script exits with an error if it is
                             absent.
    :param variable_tfs:     set of str. Gene names to force-include regardless
                             of ranking. These must already have been removed from
                             gene_df before calling this function.
    :param num_genes:        int. Maximum number of non-TF genes to include.
                             0 means include TFs only.
    :param sort_by_variance: bool. When True, rank non-TF genes by the variance
                             column (descending) before taking the top num_genes.
                             When False, use the existing row order (p-value
                             ascending).
    :return:                 list of str. Sorted (alphabetical) gene names
                             combining the selected non-TF genes and variable_tfs.
    '''
    if num_genes > 0:
        if num_genes > len(gene_df):
            # Fewer candidates than requested — include everything available.
            variable_genes_new = set(gene_df.index)
        else:
            if sort_by_variance:
                if len(gene_df.columns) < 2:
                    print("ERROR: no variance column found. Expected as the second data column "
                          "(third column of the raw CSV, after the gene-name index). Quitting")
                    sys.exit()
                var_col = gene_df.columns[1]
                gene_df = gene_df.sort_values(by=var_col, ascending=False)
            # Take the top num_genes rows after (re)sorting.
            variable_genes_new = set(gene_df.iloc[:num_genes].index.values)
    else:
        # num_genes == 0: include only the TF set, no additional non-TF genes.
        variable_genes_new = set()

    return sorted(variable_genes_new | variable_tfs)


def filter_network(
    net_df: pd.DataFrame,
    gene_set: list
) -> pd.DataFrame:
    '''
    Restrict a network edge list to edges whose both endpoints are members of
    gene_set, then remove self-loops and duplicate rows.

    The input DataFrame must have columns named "Gene1" and "Gene2" (str).
    Any additional columns (e.g. edge weights) are preserved unchanged.

    Self-loops (Gene1 == Gene2) are removed because they carry no information
    for GRN inference and are not meaningful as regulatory interactions here.

    Duplicate rows can occur in some ground-truth network sources; only the
    first occurrence of each duplicated edge is retained.

    :param net_df:   pd.DataFrame with at minimum columns "Gene1" (str) and
                     "Gene2" (str) representing directed or undirected edges.
    :param gene_set: list (or any iterable) of str gene names that define the
                     allowed node set. Edges with either endpoint absent from
                     this set are dropped.
    :return:         pd.DataFrame with the same columns as net_df, containing
                     only edges within gene_set, without self-loops or
                     duplicate rows.
    '''
    gene_set = set(gene_set)

    # Keep only edges where both endpoints are in the retained gene set.
    filtered = net_df[
        net_df.Gene1.isin(gene_set) & net_df.Gene2.isin(gene_set)
    ]

    # Remove self-loops.
    filtered = filtered[filtered.Gene1 != filtered.Gene2]

    # Remove duplicates (there are some repeated lines in the ground-truth networks!!!).
    filtered = filtered.drop_duplicates(keep='first')

    return filtered


def _process(opts: argparse.Namespace) -> None:
    '''
    Run the full input-generation pipeline for a single dataset.

    :param opts: argparse.Namespace with all required fields populated.
    '''
    include_tfs      = opts.TFs
    expr_file        = opts.expFile
    gene_ordering_file = opts.geneOrderingFile
    tf_file          = opts.TFFile
    pval_cutoff      = opts.pVal
    bf_corr          = opts.BFcorr
    num_genes        = opts.numGenes
    sort_by_variance = opts.sortByVariance

    print("\nReading %s" % expr_file)
    expr_df = pd.read_csv(expr_file, header=0, index_col=0)
    print("\nReading %s" % gene_ordering_file)
    gene_df = pd.read_csv(gene_ordering_file, header=0, index_col=0)

    if include_tfs:
        print("\nReading %s" % tf_file)
        tfs_df = pd.read_csv(tf_file, header=0)
        tfs = tfs_df[tfs_df.columns[0]]
        num_total_tfs = len(tfs)
        # limit the tfs to those present in the expression file
        tfs = tfs[tfs.isin(expr_df.index)]
        print("\t%d/%d TFs present in ExpressionData" % (len(tfs), num_total_tfs))

    # make sure the variable genes are in the Expression Data
    expr_variable_genes = set(gene_df.index.values) & set(expr_df.index.values)
    extra_genes = set(gene_df.index.values) - set(expr_df.index.values)
    if len(extra_genes) != 0:
        print("\nWARNING: %d variable genes are not in ExpressionData.csv:" % len(extra_genes))
        print(extra_genes)
        gene_df = gene_df.loc[list(expr_variable_genes)]

    # sort by the p-value column (first data column)
    pval_col = gene_df.columns[0]
    gene_df.sort_values(by=pval_col, inplace=True)

    gene_df = filter_genes_by_pvalue(gene_df, pval_cutoff, bf_corr)

    # separate TFs from the non-TF gene pool
    if include_tfs:
        tfs = tfs[tfs.isin(gene_df.index)]
        if pval_cutoff:
            print("\nIncluding %d TFs that pass the pval cutoff" % len(tfs))
        else:
            print("\nIncluding %d TFs" % len(tfs))
        variable_tfs = set(tfs)
        gene_df.drop(labels=variable_tfs, axis='index', inplace=True)
    else:
        variable_tfs = set()

    variable_genes = select_gene_set(gene_df, variable_tfs, num_genes, sort_by_variance)

    print("\nRestricting to %d genes" % len(variable_genes))
    expr_df = expr_df.loc[variable_genes]
    print("\nNew shape of Expression Data %d x %d" % (expr_df.shape[0], expr_df.shape[1]))

    out_dir = Path(opts.outPrefix)
    out_dir.mkdir(parents=True, exist_ok=True)
    expr_df.to_csv(out_dir / 'ExpressionData.csv')

    # Copy pseudotime file to output directory as PseudoTime.csv.
    pseudo_src = Path(opts.pseudoTimeFile)
    pseudo_dst = out_dir / 'PseudoTime.csv'
    pseudo_dst.write_bytes(pseudo_src.read_bytes())

    # Default stats written when no network file is provided.
    # num_tfs and num_genes fall back to the selected gene set counts.
    num_tfs_stat  = len(variable_tfs)
    num_genes_stat = len(variable_genes)
    density_stat  = float('nan')

    if opts.netFile is not None:
        net_df = pd.read_csv(opts.netFile)
        net_df = filter_network(net_df, variable_genes)
        net_df.to_csv(out_dir / 'GroundTruthNetwork.csv', index=False)
        all_nodes = set(net_df.Gene1.unique()).union(set(net_df.Gene2.unique()))
        num_tfs_stat   = expr_df[expr_df.index.isin(net_df.Gene1.unique())].shape[0]
        num_genes_stat = expr_df[expr_df.index.isin(all_nodes)].shape[0]
        density_stat   = net_df.shape[0] / ((num_tfs_stat * num_genes_stat) - num_tfs_stat)
        print("\n#TFs: %d, #Genes: %d, #Edges: %d, Density: %.3f" % (
            num_tfs_stat, num_genes_stat, net_df.shape[0], density_stat))

    # Append one row to dataset_stats.csv in the same directory as the other outputs.
    # Columns: dataset_id (outPrefix), num_tfs, num_genes, density (NaN when no network).
    stats_path = out_dir.parent / 'dataset_stats.csv'
    stats_df = pd.DataFrame([{
        'dataset_id': opts.dataset_id,
        'num_tfs':    num_tfs_stat,
        'num_genes':  num_genes_stat,
        'density':    density_stat,
    }])
    write_header = not stats_path.exists()
    stats_df.to_csv(stats_path, mode='a', header=write_header, index=False)

    print("\n\nDone with %s.\n" % opts.outPrefix)


def main():
    opts = parse_arguments()

    if opts.config is not None:
        all_opts = load_yaml_config(opts.config, opts)
    else:
        opts.outPrefix = str(Path(opts.outputDirectory) / opts.dataset_id)
        all_opts = [opts]

    for run_opts in all_opts:
        _process(run_opts)

    print("\n\nExiting...\n")


if __name__ == '__main__':
    main()
