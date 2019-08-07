import os
import math
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
import concurrent.futures
from rpy2.robjects import FloatVector
from rpy2.robjects.packages import importr
from sklearn.metrics import precision_recall_curve, roc_curve, auc

def CLR(evalObject, algorithmName):
    """
    A function to compute the infered ranked list using the
    CLR (Context Likelihood or Relatedness Network) algorithm
    from the predicted ranked edges, i.e., the outputs of
    different datasets generated from the same reference
    network, for a given algorithm.

    Parameters
    ----------
    evalObject: BLEval
      An object of class :class:`BLEval.BLEval`.

    algorithmName: str
      Name of the algorithm.


    :returns:
        - A DataFrame with containing AUPRC and AUROC values for each dataset
          for the given algorithm
    """
    evaluation_scores = []
    for dataset in tqdm(evalObject.input_settings.datasets):
        outDir = str(evalObject.output_settings.base_dir) + \
                 str(evalObject.input_settings.datadir).split("inputs")[1] + \
                 "/" + dataset["name"] + "/" + algorithmName

        rank_path = outDir + "/rankedEdges.csv"
        if not os.path.isdir(outDir):
            continue
        try:
            predDF = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
            predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
            predDF.drop_duplicates(keep = 'first', inplace=True)
            predDF.reset_index(drop = True,  inplace= True)
            predDF.EdgeWeight = predDF.EdgeWeight.abs().round(6)

            nodes = np.unique(predDF[['Gene1', 'Gene2']])
            mi_matrix = predDF.pivot(index='Gene1', columns='Gene2', values='EdgeWeight').reindex(columns=nodes, index=nodes, fill_value=0.0).values
            inferredDF = pd.DataFrame(__clr__(mi_matrix), columns=nodes, index=nodes)

            inferredDF = inferredDF.rename_axis('Gene1').reset_index().melt('Gene1', value_name='EdgeWeight', var_name='Gene2')\
              .query('Gene1 != Gene2')\
              .reset_index(drop=True)

            inferredDF.drop_duplicates(keep = 'first', inplace=True)
            inferredDF.reset_index(drop = True,  inplace= True)
            inferredDF.EdgeWeight = predDF.EdgeWeight.round(6)
            inferredDF.sort_values(by='EdgeWeight', ascending=False, inplace=False)
            inferredDF.to_csv(outDir + "/CLR-rankedEdges.csv", index=False)
            inferredDF['edge'] = inferredDF.apply(lambda x: '%s-%s' % (x.Gene1, x.Gene2), axis=1)

            # Evaluate the inferred edges - compute area under the PR and ROC curves

            ## Compute reference network
            inDir = str(evalObject.input_settings.datadir) + "/" + dataset["name"]
            refNetwork_path = inDir + "/" + dataset["trueEdges"]
            refNetwork = pd.read_csv(refNetwork_path, header=0, index_col=None)
            refNetwork['edge'] = refNetwork.apply(lambda x: '%s-%s' % (x.Gene1, x.Gene2), axis=1)
            refNetwork = refNetwork[refNetwork.Gene1!=refNetwork.Gene2]
            refNetwork['isReferenceEdge'] = 1
            refNetwork = refNetwork[["edge", "isReferenceEdge"]]

            ## Assign labels to edges in inferred network -- 1 for edges in reference
            ## network and 0 for edges not in reference network
            inferredDF = pd.merge(inferredDF, refNetwork, how='left', on='edge').fillna(0)
            inferredDF.isReferenceEdge = inferredDF.isReferenceEdge.astype(int)

            y_true = inferredDF['isReferenceEdge'].values
            y_scores = inferredDF.EdgeWeight.abs().values

            prroc = importr('PRROC')
            prCurve = prroc.pr_curve(scores_class0 = FloatVector(list(y_scores)),
                      weights_class0 = FloatVector(list(y_true)))

            prec, recall, thresholds = precision_recall_curve(y_true=y_true, probas_pred=y_scores, pos_label=1)
            fpr, tpr, thresholds = roc_curve(y_true=y_true, y_score=y_scores, pos_label=1)

            evaluation_scores.append([dataset["name"], algorithmName, prCurve[2][0], auc(fpr, tpr)])
        except:
            print("\nSkipping CLR computation for ", algorithmName, "on path", outDir)
            continue

    return pd.DataFrame(evaluation_scores, columns=["dataset", "algorithm", "AUPRC", "AUROC"])



def __clr__(mi_matrix):
    """
    A function to compute the infered network using the
    CLR (Context Likelihood or Relatedness Network) algorithm
    from a mutual information matrix.

    Parameters
    ----------
    mi_matrix: ndarray
      A 2D array storing mutual information between two nodes.

    :returns:
        - res: Infered network (a weighted adjacency matrix)
    """
    n = len(mi_matrix)
    res = np.zeros((n, n))
    avg = np.zeros(n)
    var = np.zeros(n)
    for i in range(n):
        for j in range(n):
            avg[i] += mi_matrix[i][j]
        avg[i] /= n
        for j in range(n):
            var[i] += math.pow(mi_matrix[i][j]-avg[i], 2)
        var[i] /= n

    for i in range(1, n):
        for j in range(n):
            if j < i:
                tmp = (mi_matrix[i][j] - avg[i])
                zi = 0 if (tmp < 0) else tmp*tmp/var[i]
                tmp = (mi_matrix[i][j] - avg[j])
                zj = 0 if (tmp < 0) else tmp*tmp/var[j]
                res[i][j] = math.sqrt(zi*zi+zj*zj)
                res[j][i] = math.sqrt(zi*zi+zj*zj)
    return res
