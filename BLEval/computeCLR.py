import os
import math
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
import concurrent.futures

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
        - edgelist: Infered ranked list from CLR network
    """
    for dataset in tqdm(evalObject.input_settings.datasets):
        outDir = str(evalObject.output_settings.base_dir) + \
                 str(evalObject.input_settings.datadir).split("inputs")[1] + \
                 "/" + dataset["name"] + "/" + algorithmName

        rank_path = outDir + "/rankedEdges.csv"
        if not os.path.isdir(outDir):
            continue
        try:
            predDF = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
        except:
            print("\nSkipping CLR computation for ", algorithmName, "on path", outDir)
            continue
            
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
