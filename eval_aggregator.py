#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import argparse
import eval as ev
import numpy as np
import pandas as pd
from tqdm import tqdm
import networkx as nx
from pathlib import Path
import src.plotCurves as pc
import src.computeMotifs as cm
import src.analyzeTopk as tk
from scipy.stats import spearmanr
from itertools import permutations
from collections import defaultdict
from networkx.convert_matrix import from_pandas_adjacency


import multiprocessing
from multiprocessing import Pool, cpu_count
import concurrent.futures

def computePairwiseJacc(inDict):
    jaccDF = {key:{key1:{} for key1 in inDict.keys()} for key in inDict.keys()}
    for key_i in inDict.keys():
        for key_j in inDict.keys():
            num = len(inDict[key_i].intersection(inDict[key_j]))
            den = len(inDict[key_i].union(inDict[key_j]))
            if den != 0:
                jaccDF[key_i][key_j] = num/den
            else:
                jaccDF[key_i][key_j] = 0
    return pd.DataFrame(jaccDF)

class EvalAggregator(object):
    '''
    The Evaluation object is created by parsing a user-provided configuration
    file. Its methods provide for further processing its inputs into
    a series of jobs to be run, as well as running these jobs.
    '''

    def __init__(self,
            input_settings: ev.InputSettings,
            output_settings: ev.OutputSettings) -> None:

        self.input_settings = input_settings
        self.output_settings = output_settings



    def findSpearman(self, algo_name):
        rankDict = {}
        sim_names = []
        for dataset in tqdm(self.input_settings.datasets):
            trueEdgesDF = pd.read_csv(str(self.input_settings.datadir)+'/'+ \
                                      dataset['trueEdges'],
                                sep = ',', header = 0, index_col = None)
            possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                         r = 2))
            PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}

            outDir = str(self.output_settings.base_dir) + \
                     str(self.input_settings.datadir).split("inputs")[1] + \
                     "/" + dataset["name"] + "/" + algo_name
            #algos = self.input_settings.algorithms
            rank_path = outDir+"/rankedEdges.csv"
            if not os.path.isdir(outDir):
                continue
            try:
                predEdgeDF = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
            except:
                print("Skipping spearman computation for ", algo_name, "on path", outDir)
                continue

            for key in PredEdgeDict.keys():
                subDF = predEdgeDF.loc[(predEdgeDF['Gene1'] == key.split('|')[0]) &
                               (predEdgeDF['Gene2'] == key.split('|')[1])]
                if len(subDF)>0:
                    PredEdgeDict[key] = np.abs(subDF.EdgeWeight.values[0])
            rankDict[dataset["name"]] = PredEdgeDict
            sim_names.append(dataset["name"])
        df2 = pd.DataFrame.from_dict(rankDict)
        # print(algo_name+":"+"\n\n")
        spearmanDF = df2.corr(method='spearman')
        
        
        df = spearmanDF.where(np.triu(np.ones(spearmanDF.shape),  k = 1).astype(np.bool))
    
        df = df.stack().reset_index()
        df.columns = ['Row','Column','Value']
        
        return(df.Value.median(),df.Value.mad())


    
    
    def computeJaccard(self, algo_name):
        rankDict = {}
        sim_names = []
        for dataset in tqdm(self.input_settings.datasets):
            trueEdgesDF = pd.read_csv(str(self.input_settings.datadir)+'/'+ \
                                      dataset['trueEdges'], sep = ',',
                                      header = 0, index_col = None)
            
            possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                         r = 2))
            
            TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
            PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}

            # Compute TrueEdgeDict Dictionary
            # 1 if edge is present in the ground-truth
            # 0 if edge is not present in the ground-truth
            numEdges = 0
            for key in TrueEdgeDict.keys():
                if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                       (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
                        TrueEdgeDict[key] = 1
                        numEdges += 1
                        
            outDir = str(self.output_settings.base_dir) + \
                     str(self.input_settings.datadir).split("inputs")[1] + \
                     "/" + dataset["name"] + "/" + algo_name
            
            #algos = self.input_settings.algorithms
            rank_path = outDir + "/rankedEdges.csv"
            if not os.path.isdir(outDir):
                continue
            try:
                predDF = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
            except:
                print("Skipping Jaccard computation for ", algo_name, "on path", outDir)
                continue

            predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
            predDF.drop_duplicates(keep = 'first', inplace=True)
            predDF.reset_index(drop = True,  inplace= True)
            # check if ranked edges list is empty
            # if so, it is just set to an empty set
            
            if not predDF.shape[0] == 0:

                # we want to ensure that we do not include
                # edges without any edge weight
                # so check if the non-zero minimum is
                # greater than the edge weight of the top-kth
                # node, else use the non-zero minimum value.
                predDF.EdgeWeight = predDF.EdgeWeight.round(6)
                predDF.EdgeWeight = predDF.EdgeWeight.abs()

                # Use num True edges or the number of
                # edges in the dataframe, which ever is lower
                maxk = min(predDF.shape[0], numEdges)
                edgeWeightTopk = predDF.iloc[maxk-1].EdgeWeight
                
                nonZeroMin = np.nanmin(predDF.EdgeWeight.replace(0, np.nan).values)
                bestVal = max(nonZeroMin, edgeWeightTopk)
                
                newDF = predDF.loc[(predDF['EdgeWeight'] >= bestVal)]
                rankDict[dataset["name"]] = set(newDF['Gene1'] + "|" + newDF['Gene2'])
            else:
                rankDict[dataset["name"]] = set([])
            
        Jdf = computePairwiseJacc(rankDict)
        df = Jdf.where(np.triu(np.ones(Jdf.shape),  k = 1).astype(np.bool))
        df = df.stack().reset_index()
        df.columns = ['Row','Column','Value']
        return(df.Value.median(),df.Value.mad())

    def computeEarlyPrecRec(self, algo_name):
        rankDict = {}
        sim_names = []
        for dataset in tqdm(self.input_settings.datasets):
            trueEdgesDF = pd.read_csv(str(self.input_settings.datadir)+'/'+ \
                                      dataset['trueEdges'], sep = ',',
                                      header = 0, index_col = None)
            
            possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                         r = 2))
            
            TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
            PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}

            # Compute TrueEdgeDict Dictionary
            # 1 if edge is present in the ground-truth
            # 0 if edge is not present in the ground-truth
            numEdges = 0
            trueEdges = set()
            for key in TrueEdgeDict.keys():
                if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                       (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
                        TrueEdgeDict[key] = 1
                        trueEdges.add(key)
                        numEdges += 1
                        
            outDir = str(self.output_settings.base_dir) + \
                     str(self.input_settings.datadir).split("inputs")[1] + \
                     "/" + dataset["name"] + "/" + algo_name
            
            #algos = self.input_settings.algorithms
            rank_path = outDir + "/rankedEdges.csv"
            if not os.path.isdir(outDir):
                continue
            try:
                predDF = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
            except:
                print("Skipping Jaccard computation for ", algo_name, "on path", outDir)
                continue

            predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
            predDF.drop_duplicates(keep = 'first', inplace=True)
            predDF.reset_index(drop = True,  inplace= True)
            # check if ranked edges list is empty
            # if so, it is just set to an empty set
            
            if not predDF.shape[0] == 0:

                # we want to ensure that we do not include
                # edges without any edge weight
                # so check if the non-zero minimum is
                # greater than the edge weight of the top-kth
                # node, else use the non-zero minimum value.
                predDF.EdgeWeight = predDF.EdgeWeight.round(6)
                predDF.EdgeWeight = predDF.EdgeWeight.abs()

                # Use num True edges or the number of
                # edges in the dataframe, which ever is lower
                maxk = min(predDF.shape[0], numEdges)
                edgeWeightTopk = predDF.iloc[maxk-1].EdgeWeight
                
                nonZeroMin = np.nanmin(predDF.EdgeWeight.replace(0, np.nan).values)
                bestVal = max(nonZeroMin, edgeWeightTopk)
                
                newDF = predDF.loc[(predDF['EdgeWeight'] >= bestVal)]
                rankDict[dataset["name"]] = set(newDF['Gene1'] + "|" + newDF['Gene2'])
            else:
                print(rank_path,dataset["name"])
                rankDict[dataset["name"]] = set([])
        Eprec = {}
        Erec = {}
        for dataset in tqdm(self.input_settings.datasets):
            if len(rankDict[dataset["name"]]) != 0:
                intersectionSet = rankDict[dataset["name"]].intersection(trueEdges)
                Eprec[dataset["name"]] = len(intersectionSet)/len(rankDict[dataset["name"]])
                Erec[dataset["name"]] = len(intersectionSet)/len(trueEdges)
            else:
                Eprec[dataset["name"]] = 0
                Erec[dataset["name"]] = 0
        
        return(np.median(list(Eprec.values())),Eprec)

    
    def computeSignedEarlyPrec(self, algo_name):
        rankDict = {'+':{},'-':{}}
        sim_names = []
        for sgn in ['+','-']:
            for dataset in tqdm(self.input_settings.datasets):
                trueEdgesDF = pd.read_csv(str(self.input_settings.datadir)+'/'+ \
                                          dataset['trueEdges'], sep = ',',
                                          header = 0, index_col = None)

                possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                             r = 2))


                TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
                PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}

                # Compute TrueEdgeDict Dictionary
                # 1 if edge is present in the ground-truth
                # 0 if edge is not present in the ground-truth
                numEdges = 0
                trueEdges = set()
                toRemove = []
                for key in TrueEdgeDict.keys():
                    subDF = trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene2'] == key.split('|')[1])]

                    if subDF.shape[0] > 0:
                        if  subDF['Type'].values[0] == sgn:
                            TrueEdgeDict[key] = 1
                            trueEdges.add(key)
                            numEdges += 1
                        else:
                            toRemove.append(key)
                for key in toRemove:
                    TrueEdgeDict.pop(key, None)


                outDir = str(self.output_settings.base_dir) + \
                         str(self.input_settings.datadir).split("inputs")[1] + \
                         "/" + dataset["name"] + "/" + algo_name

                #algos = self.input_settings.algorithms
                rank_path = outDir + "/rankedEdges.csv"
                if not os.path.isdir(outDir):
                    print(outDir," not found")
                    continue
                try:
                    predDF = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
                except:
                    print("Skipping signed precision computation for ", algo_name, "on path", outDir)
                    continue

                predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
                predDF.drop_duplicates(keep = 'first', inplace=True)
                predDF.reset_index(drop = True,  inplace= True)

                # Remove incorrect sign from consideration
                for idx, row in predDF.iterrows():
                    if str(row['Gene1']) + '|' + str(row['Gene2']) not in TrueEdgeDict.keys():
                        predDF.drop(idx, axis = 'index', inplace= True)
                predDF.reset_index(drop = True,  inplace= True)
                # check if ranked edges list is empty
                # if so, it is just set to an empty set

                if not predDF.shape[0] == 0:

                    # we want to ensure that we do not include
                    # edges without any edge weight
                    # so check if the non-zero minimum is
                    # greater than the edge weight of the top-kth
                    # node, else use the non-zero minimum value.
                    predDF.EdgeWeight = predDF.EdgeWeight.round(6)
                    predDF.EdgeWeight = predDF.EdgeWeight.abs()

                    # Use num True edges or the number of
                    # edges in the dataframe, which ever is lower
                    maxk = min(predDF.shape[0], numEdges)
                    edgeWeightTopk = predDF.iloc[maxk-1].EdgeWeight

                    nonZeroMin = np.nanmin(predDF.EdgeWeight.replace(0, np.nan).values)
                    bestVal = max(nonZeroMin, edgeWeightTopk)

                    newDF = predDF.loc[(predDF['EdgeWeight'] >= bestVal)]
                    rankDict[sgn][dataset["name"]] = set(newDF['Gene1'] + "|" + newDF['Gene2'])
                else:
                    print(rank_path,dataset["name"])
                    rankDict[sgn][dataset["name"]] = set([])
        Pprec = {'+':{},'-':{}}
        
        for sgn in ['+','-']:

            TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}

            # Compute TrueEdgeDict Dictionary
            # 1 if edge is present in the ground-truth
            # 0 if edge is not present in the ground-truth
            trueEdges = set()
            for key in TrueEdgeDict.keys():
                subDF = trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                       (trueEdgesDF['Gene2'] == key.split('|')[1])]

                if subDF.shape[0] > 0:
                    if  subDF['Type'].values[0] == sgn:
                        TrueEdgeDict[key] = 1
                        trueEdges.add(key)
            print(trueEdges,sgn)
            for dataset in tqdm(self.input_settings.datasets):
                if len(rankDict[sgn][dataset["name"]]) != 0 and len(trueEdges) != 0:
                    intersectionSet = rankDict[sgn][dataset["name"]].intersection(trueEdges)
                    Pprec[sgn][dataset["name"]] = len(intersectionSet)/len(rankDict[sgn][dataset["name"]])
                else:
                    Pprec[sgn][dataset["name"]] = 0
        print(algo_name, pd.DataFrame(Pprec).median(axis='index').values[0],pd.DataFrame(Pprec).median(axis='index').values[1])
        return(pd.DataFrame(Pprec).median(axis='index').values[0],pd.DataFrame(Pprec).median(axis='index').values[1])



    def create_rank_graph(self, simulations, outDir, algo_name):
        for simulation in simulations:
            simpath = outDir+"/"+simulation+"/"+algo_name
            rank_path = simpath+"/rankedEdges.csv"

            if not os.path.isdir(simpath):
                continue
            try:
                g = nx.read_weighted_edgelist(rank_path, delimiter="\t", comments="Gene1", nodetype=str)

            except FileNotFoundError:
                print("skipping graph generation for ", algo_name, "on path", outDir)


        return g

    def find_correlations(self):

        '''
        Finds the Spearman's correlation coefficient between consecutive simulations
        of the same algorithm, with the same initial parameters
        Saves the coefficients to a csv file, to be read into a Pandas dataframe
        :return:
        None
        '''
        corrDF = {}
        corrDF['Median'] = {}
        corrDF['MAD'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        print(outDir)
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            corrDF['Median'][algo[0]],corrDF['MAD'][algo[0]] = self.findSpearman(algo[0])
            
        return pd.DataFrame(corrDF)

    def find_jaccard(self):

        '''
        Computes the Jaccard Index.
        Saves the coefficients to a csv file, to be read into a Pandas dataframe.
        :return:
        None
        '''
        JaccDF = {}
        JaccDF['Median'] = {}
        JaccDF['MAD'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        print(outDir)
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            JaccDF['Median'][algo[0]], JaccDF['MAD'][algo[0]] = self.computeJaccard(algo[0])
            
        return pd.DataFrame(JaccDF)
                 
    def find_EarlyPrecRec(self):

        '''
        Computes the Jaccard Index.
        Saves the coefficients to a csv file, to be read into a Pandas dataframe.
        :return:
        None
        '''
        ePR = {}
        ePR['Early Precision'] = {}
        Eprec = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        print(outDir)
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            ePR['Early Precision'][algo[0]], Eprec[algo[0]] = self.computeEarlyPrecRec(algo[0])
        pd.DataFrame(Eprec).to_csv(outDir+'EarlyPRFull.csv')
        return pd.DataFrame(ePR)
    
    def find_SignedPrec(self):

        '''
        Computes the Jaccard Index.
        Saves the coefficients to a csv file, to be read into a Pandas dataframe.
        :return:
        None
        '''
        ePR = {}
        ePR['Activation Precision'] = {}
        ePR['Inhibition Precision'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):

            ePR['Activation Precision'][algo[0]], ePR['Inhibition Precision'][algo[0]] = self.computeSignedEarlyPrec(algo[0])
            
        return pd.DataFrame(ePR)
    
         

    def find_curves(self):

        '''
        Plot PR and ROC curves for each dataset
        for all the algorithms
        '''
        AUPRCDict = {}
        AUROCDict = {}
        uAUPRCDict = {}
        uAUROCDict = {}
        
        for dataset in tqdm(self.input_settings.datasets, 
                            total = len(self.input_settings.datasets), unit = " Dataset"):
            
            AUPRC, AUROC = pc.PRROC(dataset, self.input_settings, directed = True, selfEdges = False)
            uAUPRC, uAUROC = pc.PRROC(dataset, self.input_settings, directed = False, selfEdges = False)
            
            AUPRCDict[dataset['name']] = AUPRC
            AUROCDict[dataset['name']] = AUROC
            uAUPRCDict[dataset['name']] = uAUPRC
            uAUROCDict[dataset['name']] = uAUROC
            
        return AUPRCDict, AUROCDict, uAUPRCDict, uAUROCDict
    
    def find_motifs(self):

        '''
        compute motifs
        '''
        
        for dataset in tqdm(self.input_settings.datasets, 
                            total = len(self.input_settings.datasets), unit = " Dataset"):
            cm.motifAnalysis(dataset, self.input_settings)
    
    def analyzeTopK(self):

        '''
        compute motifs
        '''
        
        for dataset in tqdm(self.input_settings.datasets, 
                            total = len(self.input_settings.datasets), unit = " Dataset"):
            tk.outputAnalysis(dataset, self.input_settings)
    
    def time_runners(self):
        '''
        :return:
        dict of times for all dataset, for all algos
        '''
        TimeDict = dict()

        for dataset in self.input_settings.datasets:
            timevals  = self.get_time(dataset)
            TimeDict[dataset["name"]] = timevals

        return TimeDict

    def get_time(self, dataset):
        '''
        :param dataset:
        :return: dict of algorithms, along with the corresponding time for the dataset
        '''
        #TODO:modify this to incorporate multiple simulations
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/" + dataset["name"] + "/"
        algo_dict = dict()
        algos = self.input_settings.algorithms
        for algo in algos:
            path = outDir+algo[0]+"/time.txt"
            if Path(path).exists():
                time = self.parse_time_files(path)
            else:
                PTFile = str(self.input_settings.datadir.joinpath(dataset['name']+'/'+dataset['cellData']))
                PTData = pd.read_csv(PTFile,
                             header = 0, index_col = 0)
                colNames = PTData.columns
                for idx in range(len(colNames)):
                    path = outDir+algo[0]+"/time"+str(idx)+".txt"
                    time = self.parse_time_files(path)

            if time == -1:
                print("skipping time computation for ", algo[0], "on dataset", dataset["name"])
                continue

            algo_dict[algo[0]] = time

        return algo_dict

    def parse_time_files(self, path):
        """
       gets the user time given by the time command, which is stored along with the output when the algorithm runs,
       in a file called time.txt
       :return:
       float containing the time this object took to run on the dataset
        """
        try:
            with open(path, "r") as f:
                lines = f.readlines()
                line = lines[1]
                time_val = float(line.split()[-1])

        except FileNotFoundError:
            print("time output " +path+" file not found, setting time value to -1")
            time_val = -1
        except ValueError:
            print("Algorithm running failed, setting time value to -1")
            time_val = -1

        return time_val


def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('--config', default='config.yaml',
        help='Configuration file')

    return parser

def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts

def main():
    #TODO:make pretty, add comments
    opts = parse_arguments()
    config_file = opts.config

    evaluation = None

    with open(config_file, 'r') as conf:
        evaluation = ev.ConfigParser.parse(conf)
    print(evaluation)
    print('Post-Run Evaluation Summary started')

    eval_summ = EvalAggregator(evaluation.input_settings, evaluation.output_settings)
    outDir = str(eval_summ.output_settings.base_dir) + \
            str(eval_summ.input_settings.datadir).split("inputs")[1] + "/"+\
            str(eval_summ.output_settings.output_prefix) + "-"
    #TimeDict = eval_summ.time_runners()
    #pd.DataFrame(TimeDict).to_csv(outDir+'Timescores.csv')
    #AUPRCDict, AUROCDict, uAUPRCDict, uAUROCDict = eval_summ.find_curves()

    #pd.DataFrame(AUPRCDict).to_csv(outDir+'AUPRCscores.csv')
    #pd.DataFrame(AUROCDict).to_csv(outDir+'AUROCscores.csv')
    #pd.DataFrame(uAUPRCDict).to_csv(outDir+'uAUPRCscores.csv')
    #pd.DataFrame(uAUROCDict).to_csv(outDir+'uAUROCscores.csv')
    #sys.exit()
    #sPrDF = eval_summ.find_SignedPrec()
    #sPrDF.to_csv(outDir + "SingedPr.csv")
    
     # Compute correlations
    ePRDF = eval_summ.find_EarlyPrecRec()
    ePRDF.to_csv(outDir + "EarlyPR.csv")
    sys.exit()
    # Compute correlations
    jaccDict = eval_summ.find_jaccard()
    jaccDict.to_csv(outDir + "jaccData.csv")

    #corrDict = eval_summ.find_correlations()
    #corrDict.to_csv(outDir + "corrData.csv")
    #sys.exit()
    # Compute performance
    #eval_summ.analyzeTopK()





    print('Evaluation complete')


if __name__ == '__main__':
  main()
