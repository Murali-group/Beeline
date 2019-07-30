import os
import math
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
from tqdm import tqdm
from pathlib import Path
import concurrent.futures
import matplotlib.pyplot as plt
from itertools import permutations
from sklearn import preprocessing
from sklearn.metrics import precision_recall_curve, roc_curve, auc
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")


def kPRROC(evalObject, k = 1, directed = True, userReferenceNetworkFile = None):
    '''
    Computes areas under the k-precision-recall (PR) and
    and k-ROC plots for each algorithm-dataset combination.

    Parameters
    ----------
    k: int
        A predicted edge (a,b) is considered a TP if there is a
        path of length less than equal to k between a and b.
        So a 1-PR curve is just the regular PR curve.

    directed: bool
        A flag to specifiy whether to treat predictions
        as directed edges (directed = True) or
        undirected edges (directed = False).

    userReferenceNetworkFile: str
        The path to a file that specifiy reference network to be used for
        AUC calculations. Default is None. If the value is not None, the
        function will overide the reference network with the given file.
        The file should be comma separated and have following node column
        names: `Gene1` and `Gene2`.

    :returns:
        - kprrocDF: A dataframe containing k-AUPRC and k-AUROC values for each algorithm-dataset combination
    '''

    evaluation_scores = []
    for dataset in tqdm(evalObject.input_settings.datasets):
        edges = []
        for algorithmName, _ in tqdm(evalObject.input_settings.algorithms):
            outDir = str(evalObject.output_settings.base_dir) + \
                     str(evalObject.input_settings.datadir).split("inputs")[1] + \
                     "/" + dataset["name"]
            inDir = str(evalObject.input_settings.datadir) + "/" + dataset["name"]
            rank_path = outDir + "/" + algorithmName + "/rankedEdges.csv"
            if userReferenceNetworkFile is None:
                refNetwork_path = inDir + "/" + dataset["trueEdges"]
            else:
                refNetwork_path = userReferenceNetworkFile

            if not os.path.isdir(outDir) or not os.path.isdir(inDir):
                continue
            try:
                refNetwork = pd.read_csv(refNetwork_path, header=0, index_col=None)
                refNetwork = refNetwork[refNetwork.Gene1!=refNetwork.Gene2]

                all_edges_df = pd.DataFrame(list(permutations(np.unique(refNetwork.loc[:,['Gene1','Gene2']]), r = 2)), columns=['Gene1','Gene2'])

                df = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
                predEdges = pd.merge(all_edges_df, df, on=['Gene1','Gene2'], how='left')
                predEdges.EdgeWeight = predEdges.EdgeWeight.fillna(0)
                predEdges['EdgeWeight'] = predEdges['EdgeWeight'].abs()

                prec, recall, fpr, tpr, auprc, auroc = computeScores(refNetwork, predEdges, k=k, directed=directed)

                evaluation_scores.append([dataset["name"], algorithmName, auprc, auroc])
            except Exception as e:
                print(e)
                print("\nSkipping kAUC computation for %s algorithm and %s dataset on path %s" % (algorithmName, dataset["name"], outDir))
                continue

    return pd.DataFrame(evaluation_scores, columns=["dataset", "algorithm", "AUPRC", "AUROC"])


def computeScores(refNetwork, predEdgeDF, k=1):
    '''
    Computes k-precision-recall and k-ROC curves
    using scikit-learn for a given set of predictions in the
    form of a DataFrame.

    Parameters
    ----------
    refNetworkDF: DataFrame
        A pandas dataframe containing edges in the reference network.

    predEdgeDF: DataFrame
        A pandas dataframe containing the edge ranks from the prediced
        network. The indices of this dataframe are all possible edges.
        This dataframe only has one column to indicate the edge weights
        in the predicted network. Higher the weight, higher the
        edge confidence.

    k: int
        A predicted edge (a,b) is considered a TP if there is a
        path of length less than equal to k between a and b.
        So a 1-PR curve is just the regular PR curve.

    directed: bool
        A flag to specifiy whether to treat predictions
        as directed edges (directed = True) or
        undirected edges (directed = False).

    :returns:
        - prec: A list of precision values (for PR plot)
        - recall: A list of precision values (for PR plot)
        - fpr: A list of false positive rates (for ROC plot)
        - tpr: A list of true positive rates (for ROC plot)
        - AUPRC: Area under the precision-recall curve
        - AUROC: Area under the ROC curve
    '''
    G = nx.from_pandas_edgelist(refNetwork, 'Gene1', 'Gene2', True, create_using=nx.DiGraph if directed else nx.Graph)
    length = dict(nx.all_pairs_shortest_path_length(G, cutoff=k))
    if directed:
        predEdgeDF['isTP'] = predEdgeDF.apply(lambda x: 1 if (x.Gene1 in length and x.Gene2 in length[x.Gene1]) else 0, axis=1)
    else:
        predEdgeDF['isTP'] = predEdgeDF.apply(lambda x: 1 if ((x.Gene1 in length and x.Gene2 in length[x.Gene1]) or (x.Gene2 in length and x.Gene1 in length[x.Gene2])) else 0, axis=1)

    predEdgeDF.isTP = predEdgeDF.isTP.astype(int)

    y_true = predEdgeDF['isTP'].values
    y_scores = predEdgeDF['EdgeWeight'].values

    prec, recall, thresholds = precision_recall_curve(y_true=y_true, probas_pred=y_scores, pos_label=1)
    fpr, tpr, thresholds = roc_curve(y_true=y_true, y_score=y_scores, pos_label=1)

    return prec, recall, fpr, tpr, auc(recall, prec), auc(fpr, tpr)
