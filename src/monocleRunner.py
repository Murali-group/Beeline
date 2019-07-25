import os
import pandas as pd
import shutil
from pathlib import Path
import numpy as np
import sys
import subprocess


def generatePseudoTime(
        inputDir, exprDataFile="ExpressionData.csv",
        geneDataFile="GeneData.csv", cellDataFile="CellData.csv",
        out_file="PseudoTime.csv"):
    '''
    Function to generate both the PseudoTime and quantify which genes expression varies the most across time.
    '''
    outDir = inputDir.joinpath("MONOCLE")
    expr_file = inputDir.joinpath(exprDataFile)
    gene_file = outDir.joinpath(geneDataFile)
    cell_file = outDir.joinpath(cellDataFile)
    pt_file = inputDir.joinpath(out_file)
    out_file = outDir.joinpath(cellDataFile)
    time_file = outDir.joinpath("time.txt")
    processInputs(inputDir, outDir, expr_file, gene_file, cell_file)

    # TODO allow these options to be specified by the user
    lowerDetectionLimit = 0
    expressionFamily = "negbinomial.size"
    log = False

    # now run monocle
    cmdToRun = ' '.join([
        'docker run --rm -v', str(Path.cwd())+':/data/ monocle:base ',
        '/bin/sh -c \"time -v -o', "data/" + str(time_file),
        'Rscript runMonocle.R',
        '-e data/'+str(expr_file),
        '-c data/'+str(cell_file),
        '-g data/'+str(gene_file),
        '-o data/'+str(outDir),
        '-l', str(lowerDetectionLimit), 
        '-x', expressionFamily])
    if log is True:
        cmdToRun += ' --log'
    cmdToRun += '\"'
    print(cmdToRun)
    subprocess.check_call(cmdToRun, shell=True)

    # move the pseudotime file to the inputDir
    print("moving %s to %s" % (out_file, pt_file))
    shutil.move(out_file, pt_file)


def processInputs(
        inputDir, outDir, expr_file,
        gene_file, cell_file):
    '''
    Function to preprocess the expression file and generate the other necessary input files for MONOCLE
    If there are rows (genes) or columns (cells) with all 0s, they will be removed or Monocle will give an error.
    '''

    if not outDir.exists():
        print("Input folder %s does not exist, creating input folder..." % (outDir))
        outDir.mkdir(exist_ok = True)

    print("Reading expression file %s" % (expr_file))
    ExpressionData = pd.read_csv(expr_file,
                                 header = 0, index_col = 0)

    edf = ExpressionData
    rewrite = False
    # make sure there are no columns with all 0's and rows with all 0's as monocle will give an error.
    cell_exp_counts = edf.astype(bool).sum(axis=0)
    cells_with_no_exp = cell_exp_counts[cell_exp_counts == 0]
    gene_cell_counts = edf.astype(bool).sum(axis=1)
    genes_with_no_exp = gene_cell_counts[gene_cell_counts == 0]
    if len(cells_with_no_exp) != 0 or len(genes_with_no_exp) != 0:
        rewrite = True
        print("\tdropping %d cells with 0 gene counts" % (len(cells_with_no_exp)))
        print("\tdropping %d genes with 0 cell counts" % (len(genes_with_no_exp)))
        # remove those columns from the dataframe
        edf.drop(cells_with_no_exp.index, axis=1, inplace=True)
        # remove those rows from the dataframe
        edf.drop(genes_with_no_exp.index, axis=0, inplace=True)

    # make sure there aren't any duplicate gene names, as that will also cause an error
    duplicated = edf.index.duplicated(keep="first")
    duplicate_genes = edf.index[duplicated].values
    if len(duplicate_genes) > 0: 
        print("\tdropping %d rows with duplicate gene names: %s" % (len(duplicate_genes), duplicate_genes))
        print("Keeping the first instance")
        edf = edf[~duplicated]
    # do the same for the cells
    duplicated = edf.columns.duplicated(keep="first")
    duplicate_cells = edf.columns[duplicated].values
    if len(duplicate_cells) > 0: 
        print("\tdropping %d columns with duplicate cell names: %s" % (len(duplicate_cells), duplicate_cells))
        print("Keeping the first instance")
        edf = edf.loc[:,~duplicated]
    if len(duplicate_genes) > 0 or len(duplicate_cells) > 0: 
        rewrite = True

    # write the updated dataframe
    if rewrite is True:
        # append "-orig" to the original file 
        file_end = str(expr_file)[-4:]
        orig_file = str(expr_file).replace(file_end, '-orig'+file_end)
        print("moving %s to %s" % (expr_file, orig_file))
        shutil.move(expr_file, orig_file)
        # and write the new file to the same place to be used downstream
        print("writing %s" % (expr_file))
        edf.to_csv(expr_file, sep=',', header=True)

    # monocle expects a GeneData file and a CellData file. Make those here
    # this column is required by monocle
    geneDict = {}
    geneDict['gene_short_name'] = [gene.replace('x_', '') for gene in edf.index]
    geneDF = pd.DataFrame(geneDict, index = edf.index)
    print("writing %s" % (gene_file))
    geneDF.to_csv(gene_file, sep=',', header=True)

    # monocle expects a column labeled Time.
    cellDF = pd.DataFrame(pd.Series({col: 1 for col in edf.columns}), columns=["Time"])
    print("writing %s" % (cell_file))
    cellDF.to_csv(cell_file, sep=',', header=True)
