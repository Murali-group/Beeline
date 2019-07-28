import os
import math
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
import concurrent.futures
from itertools import permutations
from sklearn import preprocessing
from sklearn.metrics import precision_recall_curve, roc_curve, auc



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
                continue
        rank_df = pd.pivot_table(pd.concat(edges), values='normEdgeWeight', index='edge', columns='algo')
        rank_df['Borda'] = __normalize__(rank_df.rank(ascending=True, method=aggregationMethod).mean(axis=1).values)
        rank_df['mBorda'] = __normalize__(rank_df.rank(ascending=False, method=aggregationMethod).apply(lambda x: 1.0/(x*x)).mean(axis=1).values)
        rank_df['sBorda'] = __normalize__(rank_df[selectedAlgorithms].rank(ascending=True, method=aggregationMethod).mean(axis=1).values)
        rank_df['smBorda'] = __normalize__(rank_df[selectedAlgorithms].rank(ascending=False, method=aggregationMethod).apply(lambda x: 1.0/(x*x)).mean(axis=1).values)
        rank_df.to_csv("%s/%s-Borda.csv" % (outDir, dataset["name"]), index=False)

        refNetwork = refNetwork[["edge", "isReferenceEdge"]].set_index("edge")
        rank_df = pd.merge(rank_df, refNetwork, how='left', on='edge').fillna(0)
        rank_df.isReferenceEdge = rank_df.isReferenceEdge.astype(int)

        evaluation_scores = []
        for algo in rank_df.drop(['isReferenceEdge'], 1).columns:
            y_true = rank_df['isReferenceEdge'].values
            y_scores = rank_df[algo].values

            prec, recall, thresholds = precision_recall_curve(y_true=y_true, probas_pred=y_scores, pos_label=1)
            fpr, tpr, thresholds = roc_curve(y_true=y_true, y_score=y_scores, pos_label=1)
            evaluation_scores.append([dataset["name"], algo, auc(recall, prec), auc(fpr, tpr)])

        evaluationDFs.append(pd.DataFrame(evaluation_scores, columns=["dataset", "algorithm", "AUPRC", "AUROC"]))

    outDir = str(evalObject.output_settings.base_dir) + \
             str(evalObject.input_settings.datadir).split("inputs")[1]

    auprcDF = pd.pivot_table(pd.concat(evaluationDFs), values='AUPRC', index='algorithm', columns='dataset')
    auprcDF.to_csv("%s/%s-AUPRC-with-Borda.csv" % (outDir, dataset["name"].split("_")[0]), index=False)

    aurocDF = pd.pivot_table(pd.concat(evaluationDFs), values='AUROC', index='algorithm', columns='dataset')
    aurocDF.to_csv("%s/%s-AUROC-with-Borda.csv" % (outDir, dataset["name"].split("_")[0]), index=False)


def __normalize__(arr):
    return 1.0*(arr-np.min(arr))/(np.max(arr)-np.min(arr))
