import os
import math
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from pathlib import Path
import concurrent.futures
import matplotlib.pyplot as plt
from itertools import permutations
from sklearn import preprocessing
from sklearn.metrics import precision_recall_curve, roc_curve, auc
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")


def Borda(evalObject, selectedAlgorithms=None, aggregationMethod="average"):
    """
    A function to compute edge ranked list using the Borda method from
    the predicted ranked edges, i.e., the outputs of different datasets
    generated from the same reference network, for each dataset.
    Parameters
    ----------
    evalObject: BLEval
      An object of class :class:`BLEval.BLEval`.
    selectedAlgorithms: [str]
      List of algorithm names used to run borda method on selected
      algorithms. If nothing is provided, the function runs borda on
      all available algorithms.
    aggregationMethod: str
      Method used to aggregate rank in borda method. Available options are
      {‘average’, ‘min’, ‘max’, ‘first’}, default ‘average’
    :returns:
        - None
    """
    evaluationDFs = []
    for dataset in tqdm(evalObject.input_settings.datasets):
        edges = []
        for algorithmName, _ in tqdm(evalObject.input_settings.algorithms):
            outDir = str(evalObject.output_settings.base_dir) + \
                     str(evalObject.input_settings.datadir).split("inputs")[1] + \
                     "/" + dataset["name"]
            inDir = str(evalObject.input_settings.datadir) + "/" + dataset["name"]
            rank_path = outDir + "/" + algorithmName + "/rankedEdges.csv"
            refNetwork_path = inDir + "/" + dataset["trueEdges"]

            if not os.path.isdir(outDir) or not os.path.isdir(inDir):
                continue
            try:
                df = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
                refNetwork = pd.read_csv(refNetwork_path, header=0, index_col=None)
                refNetwork['edge'] = refNetwork.apply(lambda x: '%s-%s' % (x.Gene1, x.Gene2), axis=1)
                refNetwork = refNetwork[refNetwork.Gene1!=refNetwork.Gene2]
                refNetwork['isReferenceEdge'] = 1
                all_edges_df = pd.DataFrame(list(permutations(np.unique(refNetwork.loc[:,['Gene1','Gene2']]), r = 2)), columns=['Gene1','Gene2'])
                ranked_edges = pd.merge(all_edges_df, df, on=['Gene1','Gene2'], how='left')
                ranked_edges.EdgeWeight = ranked_edges.EdgeWeight.fillna(0)
                ranked_edges['absEdgeWeight'] = ranked_edges['EdgeWeight'].abs()
                ranked_edges['normEdgeWeight'] = pd.DataFrame(preprocessing.MinMaxScaler().fit_transform(ranked_edges[['absEdgeWeight']]))[0]
                ranked_edges['dataset'] = dataset["name"]
                ranked_edges['algo'] = algorithmName
                ranked_edges['edge'] = ranked_edges.apply(lambda x: '%s-%s' % (x.Gene1, x.Gene2), axis=1)
                edges.append(ranked_edges)
            except Exception as e:
                print("\nSkipping Borda computation for ", algorithmName, "on path", outDir)
                selectedAlgorithms.remove(algorithmName)
                continue
        rank_df = pd.pivot_table(pd.concat(edges), values='normEdgeWeight', index='edge', columns='algo')
        rank_df['BORDA'] = __normalize__(rank_df.rank(ascending=True, method=aggregationMethod).mean(axis=1).values)
        rank_df['mBORDA'] = __normalize__(rank_df.rank(ascending=False, method=aggregationMethod).apply(lambda x: 1.0/(x*x)).mean(axis=1).values)
        rank_df['sBORDA'] = __normalize__(rank_df[selectedAlgorithms].rank(ascending=True, method=aggregationMethod).mean(axis=1).values)
        rank_df['smBORDA'] = __normalize__(rank_df[selectedAlgorithms].rank(ascending=False, method=aggregationMethod).apply(lambda x: 1.0/(x*x)).mean(axis=1).values)
        rank_df['Gene1'], rank_df['Gene2'] = rank_df.index.str.split('-', 1).str
        rank_df[['Gene1','Gene2','BORDA','mBORDA','sBORDA','smBORDA']].to_csv(outDir+"/Borda.csv", index=False)



def __normalize__(arr):
    return 1.0*(arr-np.min(arr))/(np.max(arr)-np.min(arr))
