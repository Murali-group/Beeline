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
                continue
        rank_df = pd.pivot_table(pd.concat(edges), values='normEdgeWeight', index='edge', columns='algo')
        rank_df['BORDA'] = __normalize__(rank_df.rank(ascending=True, method=aggregationMethod).mean(axis=1).values)
        rank_df['mBORDA'] = __normalize__(rank_df.rank(ascending=False, method=aggregationMethod).apply(lambda x: 1.0/(x*x)).mean(axis=1).values)
        rank_df['sBORDA'] = __normalize__(rank_df[selectedAlgorithms].rank(ascending=True, method=aggregationMethod).mean(axis=1).values)
        rank_df['smBORDA'] = __normalize__(rank_df[selectedAlgorithms].rank(ascending=False, method=aggregationMethod).apply(lambda x: 1.0/(x*x)).mean(axis=1).values)
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

    evaluationDFs = pd.concat(evaluationDFs)
    outDir = str(evalObject.output_settings.base_dir) + \
             str(evalObject.input_settings.datadir).split("inputs")[1]

    auprcDF = pd.pivot_table(evaluationDFs, values='AUPRC', index='algorithm', columns='dataset')
    auprcDF.to_csv("%s/%s-AUPRC-with-Borda.csv" % (outDir, dataset["name"].split("_")[0]), index=False)

    aurocDF = pd.pivot_table(evaluationDFs, values='AUROC', index='algorithm', columns='dataset')
    aurocDF.to_csv("%s/%s-AUROC-with-Borda.csv" % (outDir, dataset["name"].split("_")[0]), index=False)

    plt.figure(figsize=(18,6))
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.6, rc={"lines.linewidth": 1.5})
    colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
              '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
              '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']
    g = sns.catplot(x="AUROC", y="algorithm", kind="box", orient="h", palette=sns.color_palette(colors),
                    order=evaluationDFs.groupby('algorithm').AUROC.median().sort_values().index.tolist(),
                    height=6, aspect=2.5, data=evaluationDFs)
    g.set(xlabel='AUROCscore', ylabel='Algorithm', title=dataset["name"].split("_")[0], xlim=(0, 1))
    sns.despine(top=False, right=False)
    plt.savefig("%s/%s-AUROCplot-with-Borda.pdf" % (outDir, dataset["name"].split("_")[0]))
    plt.savefig("%s/%s-AUROCplot-with-Borda.png" % (outDir, dataset["name"].split("_")[0]))
    plt.clf()

    plt.figure(figsize=(18,6))
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.6, rc={"lines.linewidth": 1.5})
    colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
              '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
              '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']
    g = sns.catplot(x="AUPRC", y="algorithm", kind="box", orient="h", palette=sns.color_palette(colors),
                    order=evaluationDFs.groupby('algorithm').AUROC.median().sort_values().index.tolist(),
                    height=6, aspect=2.5, data=evaluationDFs)
    g.set(xlabel='AUPRCscore', ylabel='Algorithm', title=dataset["name"], xlim=(0, 1))
    sns.despine(top=False, right=False)
    plt.savefig("%s/%s-AUPRCplot-with-Borda.pdf" % (outDir, dataset["name"].split("_")[0]))
    plt.savefig("%s/%s-AUPRCplot-with-Borda.png" % (outDir, dataset["name"].split("_")[0]))
    plt.clf()


def __normalize__(arr):
    return 1.0*(arr-np.min(arr))/(np.max(arr)-np.min(arr))
