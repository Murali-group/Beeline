import os
import pandas as pd
import numpy as np
import sys
import pandas as pd
import numpy as np
import BLRun.split_script as split
from pathlib import *

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for CNNC.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    ### ORIGINAL
    save_dir = RunnerObj.inputDir.joinpath("CNNC")
    if not save_dir.exists():
        print("Input folder for CNNC does not exist, creating input folder...")
        save_dir.mkdir(exist_ok = False)
        
    if not save_dir.joinpath("NEPDF_data/").exists():
        # input data
        expression_data = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0).T  
    
        expression_data.index = range(expression_data.shape[0])
        # expression_data.to_hdf("./expData.h5", key="expData")
        gene_labels = pd.DataFrame(expression_data.columns)
        gene_labels.columns = ['c1']
        gene_labels['c2'] = gene_labels['c1']
        # print(gene_labels)
        # gene_labels.to_csv("./geneLabels.txt", header=None, index=None, sep='\t')
        gene_pairs = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges), header = 0) # TODO maybe check for directed column
        gene_pairs.columns = ['source', 'target'] # don't care about directed
        gene_pairs = gene_pairs[['source', 'target']] # ^
        negative_pairs = split.generate_negative_samples_CNNC(df=gene_pairs, seed=1)
        negative_pairs['label'] = 0
        
        gene_pairs['label'] = 1
        gene_pairs.sort_values(by='source', axis = 0, inplace=True)
        negative_pairs.sort_values(by='source', axis=0, inplace=True)
        gene_pairs.reset_index(drop=True, inplace=True)
        negative_pairs.reset_index(drop=True, inplace=True)
        # reversed_pairs = gene_pairs.copy(deep=True)
        # reversed_pairs[['source', 'target']] = reversed_pairs[['target', 'source']]
        # reversed_pairs['label'] = 2
        # final_pairs = pd.concat([gene_pairs, reversed_pairs, negative_pairs], axis=0).sort_index(kind='merge')
        final_pairs = pd.concat([gene_pairs, negative_pairs], axis=0).sort_index(kind='merge')
        # print(gene_pairs)
        # print(negative_pairs)
        # print(final_pairs)
        # print(final_pairs.shape)
        # final_pairs.to_csv(RunnerObj.inputDir.joinpath("./genes.txt"), header=False, index=False, sep='\t')
        k = RunnerObj.params['k'] # k fold validation
        train, test = split.split_edge_cv(gene_pairs, k, RunnerObj.params['rngseed'])
        for i in range(k):   
            split.verify_split(gene_pairs, train[i], test[i], "edge")
            x = []
            y = []
            z = []
            for index, gene_pair in final_pairs.loc[test[i]].iterrows(): 
                # print(gene_pair)
                x_gene_name, y_gene_name, label = gene_pair[0], gene_pair[1], gene_pair[2]
                # print(x_gene_name)
                y.append(label)
                z.append(x_gene_name+'\t'+y_gene_name)
                # print(expression_data[x_gene_name])
                # quit()
                x_tf = np.log10(expression_data[x_gene_name] + 10 ** -2)
                x_gene = np.log10(expression_data[y_gene_name] + 10 ** -2)
                H_T = np.histogram2d(x_tf, x_gene, bins=32)
                H = H_T[0].T
                HT = (np.log10(H / (final_pairs.__len__() / 2) + 10 ** (-int(np.log10(final_pairs.__len__())))) + 4) / 4
                # print(final_pairs.__len__())
                # print((-int(np.log10(final_pairs.__len__()))))
                # print(HT)
                x.append(HT)
            if (len(x)>0):
                xx = np.array(x)[:, :, :, np.newaxis]
            else:
                xx = np.array(x)
            nepdf_dir = save_dir.joinpath('NEPDF_data')
            if not nepdf_dir.exists():
                nepdf_dir.mkdir(exist_ok=False)
            np.save(nepdf_dir.joinpath('Nxdata_tf' + str(i) + '.npy'), xx)
            np.save(nepdf_dir.joinpath('ydata_tf' + str(i) + '.npy'), np.array(y))
            np.save(nepdf_dir.joinpath('zdata_tf' + str(i) + '.npy'), np.array(z))
    ### END ORIGINAL
    return

def run(RunnerObj):
    ### ORIGINAL
    inputPath = 'data' + str(PurePath().joinpath(str(RunnerObj.inputDir).split(str(Path.cwd()))[1], "CNNC", "NEPDF_data"))
    print(inputPath)
    outDir = 'outputs' + os.path.sep + str(PurePath().joinpath(str(RunnerObj.inputDir).split("inputs/")[1], "CNNC"))
    os.makedirs(outDir, exist_ok = True)
    outPath = 'data' + os.path.sep + str(PurePath().joinpath(outDir, 'outFile.txt'))
    outPathNonFile = 'data' + os.path.sep + outDir
    cmdToRun = ' '.join(['docker run --name cnnc --gpus all --rm -v', str(Path.cwd()) + ':/data/ --expose=41269', '--entrypoint bash',
                         'cnnc:cuda -c \"/usr/bin/time -v -o', "data/" + str(outDir) + 'time.txt', # importantly, use bin/bash
                         'python cnncTrainModelKFold.py ' + str(RunnerObj.params['k']), inputPath + ' 2 ' + str(RunnerObj.params['epochs']), outPathNonFile, str(RunnerObj.params['dropouts']), str(RunnerObj.params['learning_rate']) + " > " + outPath, 
                         '\"'])
    print(cmdToRun)
    os.system(cmdToRun)

    ### END ORIGINAL
    # # print (str(Path.cwd()))
    # os.system(' '.join(['docker run --rm -v', str(Path.cwd()) + ':/data/ --expose=41269',
    #                     'cnnc:base /bin/sh -c \"time -v -o', "data/" + 'CNNCreprod/time.txt',
    #                     'python cnncDataReprod.py 13 NEPDF_data 2 > /data/CNNCreprod/outFile.txt', '\"']))


def parseOutput(RunnerObj):
    return