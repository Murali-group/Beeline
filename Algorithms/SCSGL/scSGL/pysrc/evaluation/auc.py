# scSGL - a python package for fene regulatory network inference using graph signal processing based
# signed graph learning
# Copyright (C) 2021 Abdullah Karaaslanli <evdilak@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from itertools import product, combinations, permutations

import numpy as np

from sklearn.metrics import average_precision_score, roc_auc_score

def unsigned(true_edges, pred_edges, type):

    unique_nodes = np.unique(true_edges.loc[:,['Gene1','Gene2']])
    if type == "undirected":
        possible_edges = list(combinations(unique_nodes, r = 2))
    elif type == "directed":
        possible_edges = list(permutations(unique_nodes, r = 2))
    elif type == "tfedges":
        possible_edges = set(product(set(true_edges.Gene1), set(unique_nodes)))

    true_edges_dict = {'|'.join(p):0 for p in possible_edges}
    pred_edges_dict = {'|'.join(p):0 for p in possible_edges}

    for edge in true_edges.itertuples():
        # Ignore self-edges
        if edge.Gene1 == edge.Gene2:
            continue

        if "|".join((edge.Gene1, edge.Gene2)) in true_edges_dict:
            true_edges_dict["|".join((edge.Gene1, edge.Gene2))] = 1
        
        if type == "undirected":
            if "|".join((edge.Gene2, edge.Gene1)) in true_edges_dict:
                true_edges_dict["|".join((edge.Gene2, edge.Gene1))] = 1

    for edge in pred_edges.itertuples():
        # Ignore self-edges
        if edge.Gene1 == edge.Gene2:
            continue

        if "|".join((edge.Gene1, edge.Gene2)) in pred_edges_dict:
            if np.abs(edge.EdgeWeight) > pred_edges_dict["|".join((edge.Gene1, edge.Gene2))]:
                pred_edges_dict["|".join((edge.Gene1, edge.Gene2))] = np.abs(edge.EdgeWeight)

    auprc = average_precision_score(list(true_edges_dict.values()), 
                                    list(pred_edges_dict.values()))

    auroc = roc_auc_score(list(true_edges_dict.values()), 
                          list(pred_edges_dict.values()))

    random_auprc = np.sum(np.array(list(true_edges_dict.values())))/len(true_edges_dict)
    auprc_ratio = auprc/random_auprc
    auroc_ratio = auroc/0.5

    return auprc, auroc, auprc_ratio, auroc_ratio

def signed(true_edges, pred_edges, type):
    is_pred_signed = np.count_nonzero(pred_edges.EdgeWeight >= 0) != len(pred_edges)

    unique_nodes = np.unique(true_edges.loc[:,['Gene1','Gene2']])
    if type == "undirected":
        possible_edges = list(combinations(unique_nodes, r = 2))
    elif type == "directed":
        possible_edges = list(permutations(unique_nodes, r = 2))
    elif type == "tfedges":
        possible_edges = set(product(set(true_edges.Gene1), set(unique_nodes)))

    auprc = {}
    auroc = {}
    auprc_ratio = {}
    auroc_ratio = {}
    for sgn in ["+", "-"]:
        true_edges_dict = {'|'.join(p):0 for p in possible_edges}
        pred_edges_dict = {'|'.join(p):0 for p in possible_edges}

        ignored_edges = set()
        for edge in true_edges.itertuples():
            if edge.Gene1 == edge.Gene2:
                continue

            if edge.Type == sgn:
                if "|".join((edge.Gene1, edge.Gene2)) in true_edges_dict:
                    true_edges_dict["|".join((edge.Gene1, edge.Gene2))] = 1
                
                if type == "undirected":
                    if "|".join((edge.Gene2, edge.Gene1)) in true_edges_dict:
                        true_edges_dict["|".join((edge.Gene2, edge.Gene1))] = 1
            else:
                ignored_edges.add("|".join((edge.Gene1, edge.Gene2)))

        for edge in pred_edges.itertuples():
            if edge.Gene1 == edge.Gene2:
                continue

            if is_pred_signed:
                edge_sign = "+" if edge.EdgeWeight >= 0 else "-"
                if edge_sign == sgn:
                    if "|".join((edge.Gene1, edge.Gene2)) in pred_edges_dict:
                        pred_edges_dict["|".join((edge.Gene1, edge.Gene2))] = np.abs(edge.EdgeWeight)
            else:
                if "|".join((edge.Gene1, edge.Gene2)) in ignored_edges:
                    continue
                pred_edges_dict["|".join((edge.Gene1, edge.Gene2))] = np.abs(edge.EdgeWeight)

            auprc[sgn] = average_precision_score(list(true_edges_dict.values()), 
                                    list(pred_edges_dict.values()))

            auroc[sgn] = roc_auc_score(list(true_edges_dict.values()), 
                                list(pred_edges_dict.values()))

            random_auprc = np.sum(np.array(list(true_edges_dict.values())))/len(true_edges_dict)
            auprc_ratio[sgn] = auprc[sgn]/random_auprc
            auroc_ratio[sgn] = auroc[sgn]/0.5

    return auprc, auroc, auprc_ratio, auroc_ratio
