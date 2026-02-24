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

from itertools import product, permutations
from operator import pos

import numpy as np

def unsigned(true_edges, pred_edges, type):
    true_edges_copy = true_edges.copy(deep=True)
    pred_edges_copy = pred_edges.copy(deep=True)

    # Drop self-edges and duplicates
    true_edges_copy = true_edges_copy.loc[(true_edges_copy['Gene1'] != true_edges_copy['Gene2'])]
    true_edges_copy.drop_duplicates(keep = 'first', inplace=True)
    true_edges_copy.reset_index(drop=True, inplace=True)

    pred_edges_copy = pred_edges_copy.loc[(pred_edges_copy['Gene1'] != pred_edges_copy['Gene2'])]
    pred_edges_copy.drop_duplicates(keep = 'first', inplace=True)
    pred_edges_copy.reset_index(drop=True, inplace=True)

    if type == "tfedges": # Consider only edges going out of TFs
        
        # Get a list of all possible TF to gene interactions 
        unique_nodes = np.unique(true_edges_copy.loc[:,['Gene1','Gene2']])
        possible_edges_all = set(product(set(true_edges_copy.Gene1), set(unique_nodes)))

        # Get a list of all possible interactions 
        possible_edges_no_self = set(permutations(unique_nodes, r = 2))
        
        # Find intersection of above lists to ignore self edges
        possible_edges = possible_edges_all.intersection(possible_edges_no_self)
        
        true_edges_dict = {'|'.join(p):0 for p in possible_edges}

        true_edges_str = true_edges_copy['Gene1'] + "|" + true_edges_copy['Gene2']
        true_edges_str = true_edges_str[true_edges_str.isin(true_edges_dict)]
        n_edges = len(true_edges_str)
    
        pred_edges_copy['Edges'] = pred_edges_copy['Gene1'] + "|" + pred_edges_copy['Gene2']
        # limit the predicted edges to the genes that are in the possibles
        pred_edges_copy = pred_edges_copy[pred_edges_copy['Edges'].isin(true_edges_dict)]
    else:
        true_edges_str = true_edges_copy['Gene1'] + "|" + true_edges_copy['Gene2']
        true_edges_str = set(true_edges_str.values)
        n_edges = len(true_edges_str)

    if not pred_edges_copy.shape[0] == 0:

        pred_edges_copy.loc[:, "EdgeWeight"] = pred_edges_copy.EdgeWeight.round(6).abs()
        pred_edges_copy.sort_values(by="EdgeWeight", ascending=False, inplace=True)

        # Use num True edges or the number of
        # edges in the dataframe, which ever is lower
        maxk = min(pred_edges_copy.shape[0], n_edges)
        edge_weight_topk = pred_edges_copy.iloc[maxk-1].EdgeWeight

        nnz_min = np.nanmin(pred_edges_copy.EdgeWeight.replace(0, np.nan).values)
        best_val = max(nnz_min, edge_weight_topk)

        newDF = pred_edges_copy.loc[(pred_edges_copy['EdgeWeight'] >= best_val)]
        rank = set(newDF['Gene1'] + "|" + newDF['Gene2'])

        intersectionSet = rank.intersection(true_edges_str)
        eprec = len(intersectionSet)/len(rank)
        erec = len(intersectionSet)/len(true_edges_str)

        random_eprec = n_edges/len(true_edges_dict)
        eprec_ratio = eprec/random_eprec
    else:
        eprec = 1.0
        erec = 1.0
        eprec_ratio = 1.0

    return eprec, erec, eprec_ratio